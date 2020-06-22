#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
n_cutoffs = 5000 # when thresholding graph
netw_layout = 'force-directed edgeAttribute=weight'

# NETWORK:
# edges ~ similarity
# size ~ nodes
# color ~ module color

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

# Functions.
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

#  Calculate module eigengenes.
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()
ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Evaluate relationships between modules.
adjm <- WGCNA::bicor(ME_data$eigengenes)
ne_adjm <- neten::neten(adjm)

# Cast to edge df.
edge_df <- reshape2::melt(ne_adjm)
colnames(edge_df) <- c("ModuleA","ModuleB","weight")
edge_df$ModuleA <- gsub("ME","M",edge_df$ModuleA)
edge_df$ModuleB <- gsub("ME","M",edge_df$ModuleB)
edge_df <- edge_df %>% filter(weight != 0)

# Annotate with node attributes.
noa <- data.table(Module = paste0("M",seq(1,n_modules)))

# Size ~ n nodes.
# Sizes of all modules.
module_sizes <- sapply(split(names(partition),partition),length)
names(module_sizes) <- paste0("M",names(module_sizes))
noa$Size <- module_sizes[noa$Module]

# Module color.
noa$Color <- module_colors[noa$Module]

#--------------------------------------------------------------------
## Create igraph graph objects.
#--------------------------------------------------------------------

# Creat igraph graph.
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = noa)

# Prune weak edges.
# Seq from min(edge.weight) to max to generate cutoffs.
n_edges <- length(E(g))
min_weight <- min(E(g)$weight)
max_weight <- max(E(g)$weight)
cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)

# Define a function that checks if graph is connnected at a 
# given edge weight threshold.
is_connected <- function(graph,threshold) {
filt_graph <- delete.edges(graph, 
		       which(E(graph)$weight <= threshold))
return(is.connected(filt_graph))
}

# Check if graph is connected or not at various thresholds.
checks <- sapply(cutoffs, function(threshold) {
		 is_connected(g,threshold)
		       })

# Limit is max(cutoff) at which the graph is still connected.
limit <- cutoffs[max(which(checks==TRUE))]
if (all(checks)) { stop("Error thesholding graph.") }

# Prune edges. NOTE: This removes all edge types.
g <- delete.edges(g, which(E(g)$weight <= limit))
n_edges_final <- length(E(g))

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------
# NOTE: underscores in attribute names are removed during import. 

suppressPackageStartupMessages({ library(RCy3) })

# Insure we are connected to Cytoscape.
cytoscapePing()

# Write graph to file as gml format.
# This is faster than sending to cytoscape with RCy3 function.
myfile <- file.path(netwdir, "network.gml")
write_graph(g, myfile, format = "gml")

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
	 NODE_SIZE = 35,
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

## VISUAL STYLE MAPPINGS:

# MAPPED PROPERTIES:
weight_range <- c(min(E(g)$weight), max(E(g)$weight))
size_range <- c(min(V(g)$size),max(V(g)$size))
edge_colors <- c(col2hex("gray"), col2hex("dark red"))

# List of mapped params.
mappings <- list(
NODE_FILL_COLOR = mapVisualProperty("node fill color","Color","p"),
NODE_LABEL = mapVisualProperty("node label","symbol", "p"),
EDGE_TRANSPARENCY = mapVisualProperty("edge transparency", "weight", 
			      "c", weight_range, c(155, 255)),
EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty("edge stroke unselected paint",
					 "weight", "c",weight_range,
					 edge_colors),
NODE_SIZE = mapVisualProperty("node size", "size", "c", size_range, c(35, 100)))

# Create a visual style.
createVisualStyle(style.name, defaults = defaults, 
		  mappings = mappings)

# Apply to graph.
setVisualStyle(style.name); Sys.sleep(3) 

# Apply layout.
layoutNetwork(netw_layout); Sys.sleep(2); fitContent()

# Free up some memory.
cytoscapeFreeMemory()

#---------------------------------------------------------------------
## Save
#---------------------------------------------------------------------

# When done, save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Module_Network.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

message("\nDone!")
