#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
n_cutoffs = 5000 # number of values to try when thresholding graph
netw_layout = 'force-directed edgeAttribute=weight'

#--------------------------------------------------------------------
## Misc functions.
#--------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Define a function that checks if graph is connnected at a 
# given edge weight threshold.
is_connected <- function(graph,threshold,ncomp=1) {
	filt_graph <- igraph::delete.edges(graph, 
				   which(igraph::E(graph)$weight <= threshold))
	return(components(filt_graph)$no == ncomp)
}

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
  library(RCy3) # For talking to Cytoscape.
  library(dplyr) # For manipulating data.
  library(igraph) # For creating graphs.
  library(data.table) # For working with tables.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
netwdir <- file.path(root, "networks")
figsdir <- file.path(root, "figs","Modules")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the data from root/data.
data(tmt_protein)

# Load the graph partition:
data(partition)

# Load module color assignments.
data(module_colors)

# Module names are "M#"
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# Cast TMT data into a matrix.
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()

# Compute eigengenes.
ME_data <- WGCNA::moduleEigengenes(dm,colors=partition)
ME_dm <- ME_data$eigengenes
colnames(ME_dm) <- gsub("ME","M",colnames(ME_dm))
rownames(ME_dm) <- gsub("ME","M",rownames(ME_dm))

# Compute bicor.
adjm <- WGCNA::bicor(ME_dm)

# SUBSET and then enhance!
ne_adjm <- neten::neten(adjm[sig_modules,sig_modules])

#--------------------------------------------------------------------
## Create igraph object.
#--------------------------------------------------------------------

# Cast adjm to edges df.
edge_df <- reshape2::melt(ne_adjm)
colnames(edge_df) <- c("NodeA","NodeB","weight")
# Remove edges of weight 0.
edge_df <- filter(edge_df,weight != 0) 

# Node attributes.
noa <- data.table(Node = sig_modules)
# Size ~ N nodes contained by modules within the community.
noa$Size <- sapply(modules[noa$Node],length)
noa$Color <- module_colors[noa$Node]
# Filter.
noa <- noa %>% filter(Node %in% sig_modules)

# Create igraph graph.
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = noa)

#--------------------------------------------------------------------
## Threshold graph.
#--------------------------------------------------------------------

ncomp = 1
# Prune weak edges.
# Seq from min(edge.weight) to max to generate cutoffs.
n_edges <- length(E(g))
min_weight <- min(E(g)$weight)
max_weight <- max(E(g)$weight)
cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)
# Check if graph is connected or not at various thresholds.
checks <- sapply(cutoffs, function(threshold) {
		 is_connected(g,threshold,ncomp)
		       })
# What to do if graph is already unconnected?
if(!any(checks)) { stop("Graph is unconnected.") }
# Limit is max(cutoff) at which the graph is still connected.
limit <- cutoffs[max(which(checks==TRUE))]
# Prune edges. NOTE: This removes all edge types.
subg <- delete.edges(g, which(E(g)$weight <= limit))
n_edges_final <- length(E(subg))
# Status.
dt <- data.table("Initial Number of Edges"=n_edges,
		 "Final Number of Edges"=n_edges_final)
knitr::kable(dt)

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------
# NOTE: underscores in attribute names are removed during import. 

# Insure we are connected to Cytoscape.
cytoscapePing()

generate_graph()

generate_graph <- function() {
# Write graph to file as gml format.
# This is faster than sending to cytoscape with RCy3 function.
myfile <- file.path(netwdir, "network.gml")
write_graph(subg, myfile, format = "gml")
# Load into Cytoscape.
winfile <- gsub("/mnt/d/", "D:/", myfile)
cysnetw <- importNetworkFromFile(winfile)
Sys.sleep(2); unlink(myfile) # clean-up
# VISUAL STYLE DEFAULTS:
# Create a visual style.
style.name = "myStyle"
# VISUAL STYLE DEFAULTS:
defaults <- list(
	 NETWORK_TITLE = "RCy3 Network",
	 NODE_FILL_COLOR = col2hex("gray"),
	 NODE_TRANSPARENCY = 200,
	 NODE_SIZE = 15,
	 NODE_SHAPE = "ellipse",
	 NODE_LABEL_TRANSPARENCY = 255,
	 NODE_LABEL_FONT_SIZE = 12,
	 NODE_LABEL_COLOR = col2hex("black"),
	 NODE_BORDER_TRANSPARENCY = 200,
	 NODE_BORDER_WIDTH = 4,
	 NODE_BORDER_PAINT = col2hex("black"),
	 NODE_TRANSPARENCY = 200,
	 EDGE_STROKE_UNSELECTED_PAINT = col2hex("black"),
	 EDGE_WIDTH = 2,
	 NETWORK_BACKGROUND_PAINT = col2hex("white")
	)
# MAPPED VISUAL PROPERTIES:
weight_range <- c(min(E(g)$weight), max(E(g)$weight))
size_range <- c(min(V(g)$Size),max(V(g)$Size))
edge_colors <- c(col2hex("gray"), col2hex("dark red"))
# List of mapped params.
mappings <- list(
NODE_FILL_COLOR = mapVisualProperty("node fill color","Color","p"),
NODE_LABEL = mapVisualProperty("node label","name", "p"),
EDGE_TRANSPARENCY = mapVisualProperty("edge transparency", "weight", 
			      "c", weight_range, c(155, 255)),
EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty("edge stroke unselected paint",
					 "weight", "c",weight_range,
					 edge_colors),
NODE_SIZE = mapVisualProperty("node size", "Size", "c", size_range, c(35, 150)))
# Create a visual style.
createVisualStyle(style.name, defaults = defaults, 
		  mappings = mappings)
# Apply to graph.
setVisualStyle(style.name); Sys.sleep(3) 
# Apply layout.
layoutNetwork(netw_layout)
Sys.sleep(2); fitContent()
# Free up some memory.
cytoscapeFreeMemory()
}

# Save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Module_Network.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

