#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
n_cutoffs = 1000 # number of values to try when thresholding graph
netw_layout = 'force-directed edgeAttribute=weight'

# Create a network showing relationships between modules.
# In this network edges will be ~ the similarity between two modules 
# summary profiles. Node size ~ the number of nodes in a module, and
# node color ~ module color

#--------------------------------------------------------------------
## Misc functions.
#--------------------------------------------------------------------

# Define a function that checks if graph is connnected at a 
# given edge weight threshold.
is_connected <- function(graph,threshold) {
	filt_graph <- delete.edges(graph, 
				   which(E(graph)$weight <= threshold))
	return(is.connected(filt_graph))
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
part <- paste0("M",partition)
names(part) <- names(partition)
partition <- part

# Cast data into a matrix.
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()

#--------------------------------------------------------------------
## Create adjacency matrix.
#--------------------------------------------------------------------

# Calculate module eigengenes -- a summary of a module.
ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Evaluate relationships between modules with bicor.
adjm <- WGCNA::bicor(ME_data$eigengenes)

# Perform network enhancement.
ne_adjm <- neten::neten(adjm)

#--------------------------------------------------------------------
## Create igraph object.
#--------------------------------------------------------------------

# Cast to edge df.
edge_df <- reshape2::melt(ne_adjm)
colnames(edge_df) <- c("ModuleA","ModuleB","weight")
edge_df$ModuleA <- gsub("ME","",edge_df$ModuleA) # fix names
edge_df$ModuleB <- gsub("ME","",edge_df$ModuleB) # fix names

# Remove edges of weight 0 and any edges to M0.
edge_df <- edge_df %>% filter(weight != 0) %>% 
	filter(ModuleA != "M0") %>% filter(ModuleB != "M0")

# Node attributes.
noa <- data.table(Module = unique(partition)[order(unique(partition))])
noa <- noa %>% filter(Module != "M0")

# Size ~ n nodes.
node_size <- sapply(split(names(partition),partition),length)
noa$Size <- node_size[noa$Module]

# Module color.
noa$Color <- module_colors[noa$Module]

# Creat igraph graph.
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = noa)

#--------------------------------------------------------------------
## Threshold graph.
#--------------------------------------------------------------------

# Prune weak edges.
# Seq from min(edge.weight) to max to generate cutoffs.
n_edges <- length(E(g))
min_weight <- min(E(g)$weight)
max_weight <- max(E(g)$weight)
cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)
message(paste("\nInitial number of edges:",formatC(n_edges,big.mark=",")))

# Check if graph is connected or not at various thresholds.
message("\nThresholding graph.")
checks <- sapply(cutoffs, function(threshold) {
		 is_connected(g,threshold)
		       })

# What to do if graph is already unconnected?
if(!any(checks)) { stop("Graph is unconnected.") }

# Limit is max(cutoff) at which the graph is still connected.
limit <- cutoffs[max(which(checks==TRUE))]

# Prune edges. NOTE: This removes all edge types.
g <- delete.edges(g, which(E(g)$weight <= limit))
n_edges_final <- length(E(g))
message(paste("\nFinal number of edges:",formatC(n_edges_final,big.mark=",")))

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------
# NOTE: underscores in attribute names are removed during import. 

# Insure we are connected to Cytoscape.
cytoscapePing()

# Write graph to file as gml format.
# This is faster than sending to cytoscape with RCy3 function.
myfile <- file.path(netwdir, "network.gml")
write_graph(g, myfile, format = "gml")

# Load into Cytoscape.
winfile <- gsub("/mnt/d/", "D:/", myfile)
cysnetw <- importNetworkFromFile(winfile); Sys.sleep(2)
unlink(myfile) # clean-up

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
layoutNetwork(netw_layout); Sys.sleep(2); fitContent()

# Free up some memory.
cytoscapeFreeMemory()

# Save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Module_Network.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

#--------------------------------------------------------------------
## Alternative approach: plot major communities seperately.
#--------------------------------------------------------------------

# Creat igraph graph.
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = noa)

# Cluster with louvain into large communities.
louvain_part <- cluster_louvain(g, weights = E(g)$weight)
lpartition <- louvain_part$membership
names(lpartition) <- louvain_part$names

message(paste("\nNumber of Louvain Communities:", length(communities)))
message(paste("\nLouvain Modularity:",round(louvain_part$modularity[1],4)))

# All louvain communities.
communities <- split(names(lpartition),lpartition)
names(communities) <- paste0("C",names(communities))

## Loop to generate community graphs.
for (i in c(1:length(communities))){

	# Subset graph.
	subg <- induced_subgraph(g, communities[[i]])
	
	# Threshold.
	message("\nThresholding graph.")
	n_edges <- length(E(subg))
	min_weight <- min(E(subg)$weight)
	max_weight <- max(E(subg)$weight)
	cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)
	message(paste("\nInitial number of edges:",formatC(n_edges,big.mark=",")))
	checks <- sapply(cutoffs, function(threshold) {
				 is_connected(subg,threshold)
		  })
	# Check graph should not be unconnected.
	if(!any(checks)) { stop("Graph is unconnected.") }
	limit <- cutoffs[max(which(checks==TRUE))]
	subg <- delete.edges(subg, which(E(subg)$weight <= limit))
	n_edges_final <- length(E(subg))
	message(paste("\nFinal number of edges:",formatC(n_edges_final,big.mark=",")))

	# Create Cytoscape graphs.
	myfile <- file.path(netwdir, paste0("Community_",i,".gml"))
	write_graph(subg, myfile, format = "gml")
	winfile <- gsub("/mnt/d/", "D:/", myfile)
	cysnetw <- importNetworkFromFile(winfile)
	Sys.sleep(2); unlink(myfile) # clean-up

	# Create a visual style.
	style.name = paste("community_",i,".style")
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
	weight_range <- c(min(E(subg)$weight), max(E(subg)$weight))
	size_range <- c(min(V(subg)$Size),max(V(subg)$Size))
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
	createVisualStyle(style.name, defaults = defaults, 
			  mappings = mappings)

	# Apply syle and layout to graph.
	setVisualStyle(style.name); Sys.sleep(3) 
	layoutNetwork(netw_layout); Sys.sleep(2); fitContent()

	# Free up some memory.
	cytoscapeFreeMemory(); Sys.sleep(2)
}

# Save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Module_Network.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)
