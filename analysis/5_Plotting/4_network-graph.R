#!/usr/bin/env Rscript

# title: SwipProteomics
# description: generate an overview graph of the network
# author: Tyler W Bradshaw

## ---- input
save_image = FALSE
root = "~/projects/SwipProteomics"


## ---- functions

mask <- function(dm, threshold) {
	# Apply a mask (threshold) to a matrix.
       	mask_ <- matrix(as.numeric(dm > threshold),nrow=nrow(dm), ncol=nrow(dm))
	return(mask_ * dm)
}


is_connected <- function(adjm,weighted=TRUE,diag=FALSE,mode="undirected",...) {
	# Check if the graph of an adjacency matrix is connected.
	g <- graph_from_adjacency_matrix(adjm,
					 weighted=weighted,
					 diag=diag,
					 mode=mode,...)
	return(is.connected(g))
}


## ---- Set-up the workspace

# Load renv
renv::load(root, quiet=TRUE)

# Global imports
suppressPackageStartupMessages({
  library(RCy3) # For talking to Cytoscape
  library(dplyr) # For manipulating data
  library(igraph) # For creating graphs
  library(data.table) # For working with tables
})

# project specific functions and data
devtools::load_all(root, quiet=TRUE)

# Project directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
figsdir <- file.path(root, "figs","Networks")

# Output directory for cytoscape networks
netwdir <- file.path(root,"networks")
if (!dir.exists(netwdir)) {
	dir.create(netwdir)
}

## ---- Load the data

data(module_colors)
data(ne_surprise_partition)

# ne_adjm in root/rdata
myfile <- file.path(root,"rdata","ne_adjm.rda")
load(myfile)


## ---- Create igraph graph

# Drop M0 from graph.
M0 <- names(which(partition==0))
idx <- colnames(ne_adjm) %in% M0
sub_adjm <- ne_adjm[!idx,!idx]

# Threshold the graph
# Search for an appropriate threshold
# By manual search the 'best' threshold is ...
#is_connected(mask(sub_adjm, threshold = 35.034))

# Create igraph graph from thresholded adjm
threshold = 35.034
stopifnot(is_connected(mask(sub_adjm, threshold)))

adjm <- mask(sub_adjm, threshold)
g <- graph_from_adjacency_matrix(adjm, mode="undirected", diag=F, weighted=T)

# Stop if graph is not connected
if (!is.connected(g)) { stop("The graph is not connected!") }

# Add module and color attributes
g <- set_vertex_attr(g,"Module",value=paste0("M",partition[names(V(g))]))
g <- set_vertex_attr(g,"Color",value=module_colors[V(g)$Module])

# Cys Defaults
netw_layout='force-directed edgeAttribute=weight' 

# Check that we are connected to Cytoscape
cytoscapePing()

# Write graph to file this is faster than sending to cytoscape
myfile <- file.path(netwdir, paste0("network", ".gml"))
write_graph(g, myfile, format = "gml")

# Send to Cytoscape
# NOTE: underscores in attribute names are removed
winfile <- gsub("/mnt/d/|\\~/", "D:/", myfile)
cysnetw <- importNetworkFromFile(winfile)

Sys.sleep(3)
unlink(myfile)

# For large networks you need to tell cytoscape to create a network view
result = tryCatch({
	getNetworkViews()
}, warning = function(w) {
	message(w)
}, error = function(e) {
	commandsPOST("view create")
}, finally = {
	message("Created Network View!")
	Sys.sleep(3)
})

# Visual Property defaults:
defaults <- list(
	 NETWORK_TITLE = "network",
	 NODE_FILL_COLOR = col2hex("gray"),
	 NODE_TRANSPARENCY = 200,
	 NODE_SIZE = 35,
	 NODE_SHAPE = "ellipse",
	 NODE_LABEL_TRANSPARENCY = 0,
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

# MAPPED PROPERTIES:
weight_range <- range(E(g)$weight)

# List of mapped params.
edge_colors = c(col2hex("gray"), col2hex("dark red"))
mappings <- list(
NODE_FILL_COLOR = mapVisualProperty("node fill color","Color","p"),
EDGE_TRANSPARENCY = mapVisualProperty("edge transparency", "weight", 
			      "c", weight_range, c(155, 255)),
EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty("edge stroke unselected paint",
					 "weight", "c",weight_range,
					 edge_colors))

# Create a visual style.
createVisualStyle("mysteez", defaults = defaults, 
		  mappings = mappings)

# Apply to graph.
invisible({ setVisualStyle("mysteez") })
Sys.sleep(3)

# Apply layout
invisible({ layoutNetwork(netw_layout) })
Sys.sleep(3)
fitContent()

# Save network image
if (save_image) {
  netw_image <- file.path(figsdir, "Network_Overview")
  winfile <- gsub("/mnt/d/", "D:/", netw_image)
  exportImage(winfile, "svg")
  message("\nConvert svg image to tiff before pushing to git (too big)!")
}

# Free up some memory
invisible({ cytoscapeFreeMemory() })

# Save cytoscape sesh.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,"Network")
winfile <- gsub("/mnt/d/|\\~/","D:/",myfile) 
saveSession(winfile)

message("\nDone!")
