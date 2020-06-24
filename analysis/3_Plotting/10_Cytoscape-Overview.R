#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Plot an overview of the network.

#----------------------------------------------------------------------
## Misc functions.
#----------------------------------------------------------------------

mask <- function(dm, threshold) {
	# Apply a threshold (mask) a matrix.
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

# Functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
figsdir <- file.path(root, "figs","Networks")

# Output directory for cytoscape networks.
netwdir <- file.path(root,"networks")
if (!dir.exists(netwdir)) {
	dir.create(netwdir)
}

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the graph partition:
data(partition)

# Load the colors.
data(module_colors)

# Load networks.
data(ne_adjm) # loads "edges", then cast to adjm.
ne_adjm <- convert_to_adjm(edges)

#--------------------------------------------------------------------
## Create igraph graph.
#--------------------------------------------------------------------

# Threshold the graph.
# Search for an appropriate threshold.
# By manual search the 'best' threshold is somewhere between 1.5 and 2.0.
#is_connected(mask(ne_adjm, threshold = 1.6051))

# Create igraph graph from thresholded adjm.
threshold = 1.6051
adjm <- mask(ne_adjm, threshold)
g <- graph_from_adjacency_matrix(adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)

# Stop if graph is not connected.
if (!is.connected(g)) { stop("The graph is not connected!") }

# Add module and color attributes.
g <- set_vertex_attr(g,"Module",value=paste0("M",partition[names(V(g))]))
g <- set_vertex_attr(g,"Color",value=module_colors[V(g)$Module])

# Cys Defaults.
netw_layout='force-directed edgeAttribute=weight' 

# Check that we are connected to Cytoscape.
cytoscapePing()

# Write graph to file this is faster than sending to cytoscape.
myfile <- file.path(netwdir, paste0("network", ".gml"))
write_graph(g, myfile, format = "gml")

# Send to Cytoscape.
# NOTE: underscores in attribute names are removed. 
winfile <- gsub("/mnt/d/", "D:/", myfile)
<<<<<<< HEAD
cysnetw <- importNetworkFromFile(winfile); Sys.sleep(2)
unlink(myfile)

# For REALLY large networks you need to tell cytoscape to create a view of the
# network. I coded this before...
result = tryCatch({
	getNetworkViews()
}, warning = function(w) {
	print(w)
}, error = function(e) {
	commandsPOST("view create")
}, finally = {
	print("Created Network View!")
	Sys.sleep(2)
})
=======
cysnetw <- importNetworkFromFile(winfile)
Sys.sleep(2); unlink(myfile)

# For REALLY large networks you need to tell cytoscape to create a view of the
# network. I coded this before...

# Check if network view was successfully created.
if (length(cysnetw)==0) {
	result = tryCatch({
		getNetworkViews()
	}, warning = function(w) {
		print(w)
	}, error = function(e) {
		commandsPOST("view create")
	}, finally = {
		print("Created Network View!")
	})
} Sys.slee(2)
>>>>>>> dc80f240c117e84f0c43892d3655ed3b3235e4ec

# Visual Property defaults:
defaults <- list(
	 NETWORK_TITLE = "THIS IS A NETWORK",
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
setVisualStyle("mysteez"); Sys.sleep(3)

# Apply layout.
layoutNetwork(netw_layout); Sys.sleep(2); fitContent()

# Save network image.
netw_image <- file.path(figsdir, "Network_Overview")
winfile <- gsub("/mnt/d/", "D:/", netw_image)
exportImage(winfile, "svg")

# Convert SVG to tiff.
<<<<<<< HEAD
#system(command = paste("svg2tiff",paste0(netw_image,".svg")))
=======
system(command = paste("svg2tiff",paste0(netw_image,".svg")))
>>>>>>> dc80f240c117e84f0c43892d3655ed3b3235e4ec

# Free up some memory.
cytoscapeFreeMemory()

<<<<<<< HEAD
# Save cys session. 
# NOTE: NETWORK IS TOO BIG!
#myfile <- file.path(netwdir,"Network_Overview.cys")
#winfile <- gsub("/mnt/d/","D:/",myfile) 
#saveSession(winfile)
=======
# Save cys session.
myfile <- file.path(netwdir,"Network_Overview.cys")
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)
>>>>>>> dc80f240c117e84f0c43892d3655ed3b3235e4ec
