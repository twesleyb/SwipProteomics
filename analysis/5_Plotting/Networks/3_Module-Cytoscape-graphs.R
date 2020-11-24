#!/usr/bin/env Rscript

# title: SwipProteomics
# description: generate cytoscape networks
# authors: Tyler W Bradshaw

## Analysis options:
file_prefix  = ""

## function: create a cytoscpae graph of each module --------------------------

createCytoscapeGraph <- function(netw_g, ppi_g, nodes, module_name,
			  n_cutoffs=5000, netwdir=getwd(),
			  imgsdir=getwd(), 
			  netw_layout='force-directed edgeAttribute=weight') {
	# Imports
	suppressPackageStartupMessages({
		library(RCy3)
	})
	# Default layout.
	if (is.null(netw_layout)) { 
		netw_layout <- 'force-directed edgeAttribute=weight' 
	}
	# Subset graph: keep nodes in module.
	graph <- netw_g
	idx <- match(nodes, names(V(graph)))
	subg <- induced_subgraph(graph, vids = V(graph)[idx])
	# Node Size ~ hubbiness or importance in its subgraph.
	adjm <- as.matrix(as_adjacency_matrix(subg,attr="weight"))
	node_importance <- apply(adjm,2,sum)
	subg <- set_vertex_attr(subg,"size",value=node_importance[names(V(subg))])
	# Prune weak edges.
	# Seq from min(edge.weight) to max to generate cutoffs.
	n_edges <- length(E(subg))
	min_weight <- min(E(subg)$weight)
	max_weight <- max(E(subg)$weight)
	n_cutoffs = 5000
	cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)
	# Define a function that checks if graph is connnected at a 
	# given edge weight threshold.
	is_connected <- function(graph,threshold) {
		filt_graph <- delete.edges(graph, 
				       which(E(graph)$weight <= threshold))
		return(is.connected(filt_graph))
	}
	# Check if graph is connected or not at various thresholds.
	# NOTE: this can take a little time if n is high
	checks <- sapply(cutoffs, function(threshold) {
				 is_connected(subg,threshold)
				       })
	# Limit is max(cutoff) at which the graph is still connected.
	limit <- cutoffs[max(which(checks==TRUE))]
	if (all(checks)) { stop("Error thesholding graph.") }
	# Prune edges. NOTE: This removes all edge types.
	g <- delete.edges(subg, which(E(subg)$weight <= limit))
	n_edges <- length(E(g))
	# Write graph to file this is faster than sending to cytoscape.
	myfile <- file.path(netwdir, paste0(module_name, ".gml"))
	write_graph(g, myfile, format = "gml")
	# Send to Cytoscape.
	# NOTE: underscores in attribute names are removed. 
	winfile <- gsub("/mnt/d/", "D:/", myfile)
	cysnetw <- importNetworkFromFile(winfile)
	# keep getting errors with keys -- make sure they are there!
	Sys.sleep(2); unlink(myfile)
	####################################################################
	## VISUAL STYLE DEFAULTS:
	####################################################################
	# Create a visual style.
	style.name <- paste(module_name, "style", sep = "-")
	# VISUAL STYLE DEFAULTS:
	defaults <- list(
		 NETWORK_TITLE = "THIS IS A NETWORK",
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
	####################################################################
	## VISUAL STYLE MAPPINGS:
	####################################################################
# MAPPED PROPERTIES:
weight_range <- c(min(E(g)$weight), max(E(g)$weight))
size_range <- c(min(V(g)$size),max(V(g)$size))
edge_colors <- c(col2hex("gray"), col2hex("dark red"))
# List of mapped params.
mappings <- list(
NODE_FILL_COLOR = mapVisualProperty("node fill color","Color","p"),
NODE_LABEL = mapVisualProperty("node label","Protein", "p"),
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
	setVisualStyle(style.name)
	Sys.sleep(3)
	# Collect PPI edges.
	idx <- match(nodes, names(V(ppi_g)))
		subg <- induced_subgraph(ppi_g, 
					 vids = V(ppi_g)[idx])
	edge_list <- apply(as_edgelist(subg, names = TRUE), 1, as.list)
	# If edge list is only of length 1, unnest it to avoid problems.
	if (length(edge_list) == 1) {
		edge_list <- unlist(edge_list, recursive = FALSE)
	}
	# Add PPI edges to Cytoscape graph.
	if (length(edge_list) > 0) {
		ppi_edges <- addCyEdges(edge_list)
		# Add PPIs and set to black.
		selected_edges <- selectEdges(ppi_edges, by.col = "SUID")
		# Set to black with edge bend.
		namen <- "EDGE_STROKE_UNSELECTED_PAINT"
		setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = col2hex("black"), 
				      visual.property = namen,
				      bypass = TRUE
				      )
		setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = TRUE,
				      visual.property = "EDGE_BEND",
				      bypass = TRUE
				      )
	} # Ends IF statement.
	#  Clean-up.
	clearSelection()
	Sys.sleep(2) 
	# Apply layout.
	layoutNetwork(netw_layout)
	Sys.sleep(2)
	fitContent()
	# Mask color of non-significant nodes.
	sig <- names(V(g))[V(g)$sigprot == 1]
	ns <- names(V(g))[names(V(g)) %notin% sig]
	if (length(ns) > 0) {
		setNodeColorBypass(ns,new.colors=col2hex("gray"))
	}
	# Bold border of BioID proteins.
	wash_nodes <- names(V(g))[V(g)$isWASH==1]
	if (length(wash_nodes) > 0) {
		setNodeBorderWidthBypass(wash_nodes,new.sizes=10)
	}
	file_prefix <- formatC(as.numeric(gsub("M","",module_name)),
			       width=3,flag=0)
	netw_image <- file.path(imgsdir, paste(file_prefix,module_name,sep="_"))
	winfile <- gsub("/mnt/d/", "D:/", netw_image)
	exportImage(winfile, "svg")
	# Free up some memory.
	cytoscapeFreeMemory()
} #EOF


## ---- Set-up the workspace

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# global imports
suppressPackageStartupMessages({
  library(RCy3) # For talking to Cytoscape
  library(dplyr) # For manipulating data
  library(igraph) # For creating graphs
  library(data.table) # For working with tables
})

# project functions and data
suppressWarnings({ devtools::load_all() })

# project directories
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

# Load the data from root/data
data(gene_map)
data(sigprots)
data(partition)
data(msstats_prot)
data(module_colors)
data(msstats_results)
data(wash_interactome); wash_prots <- wash_interactome


# Load networks
myfile <- file.path(root,"rdata","ne_adjm.rda")
load(myfile) # ne_adjm
myfile <- file.path(root,"rdata","adjm.rda")
load(myfile) # adjm
myfile <- file.path(root,"rdata","ppi_adjm.rda")
load(myfile) # ppi_adjm


## ---- Create igraph graph objects

# Create a list of all modules. 
module_list <- split(names(partition),partition)[-1] # drop M0
names(module_list) <- paste0("M",names(module_list))

# Insure that matrices are in matching order
check <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!check) { stop("input adjacency matrices should be of matching dimensions") }

# Create igraph graph objects
# NOTE: graph edge weight is enhanced(cor)
netw_g <- graph_from_adjacency_matrix(ne_adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
ppi_g <- graph_from_adjacency_matrix(ppi_adjm,mode="undirected",diag=FALSE,
				     weighted=TRUE)

# Annotate graph's with protein names.
proteins <- toupper(gene_map$symbol[match(names(V(netw_g)),gene_map$uniprot)])
netw_g <- set_vertex_attr(netw_g,"protein",value = proteins)
proteins <- toupper(gene_map$symbol[match(names(V(ppi_g)),gene_map$uniprot)])
ppi_g <- set_vertex_attr(ppi_g,"protein",value = proteins)


#--------------------------------------------------------------------
## Annotate graphs with additional meta data.
#--------------------------------------------------------------------

# Collect meta data from tmt_protein.
tmp_dt <- data.table(Protein = names(V(netw_g)),
		  Module = paste0("M",partition[names(V(netw_g))]))
noa <- left_join(tmp_dt, msstats_results, by = "Protein") %>% 
	filter(Contrast == "Mutant-Control")

# Add module colors.
noa$Color <- module_colors[noa$Module]

# Add WASH annotation.
noa$isWASH <- as.numeric(noa$Protein %in% wash_prots)

# Add sig prot annotations.
noa$sigprot <- as.numeric(noa$Protein %in% sigprots)

# Loop to add node attributes to netw_graph.
for (i in c(1:ncol(noa))) {
	namen <- colnames(noa)[i]
	col_data <- setNames(noa[[i]],nm=noa$Protein)
	netw_g <- set_vertex_attr(netw_g,namen,value=col_data[names(V(netw_g))])
}

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------

# Network images will be saved in networks/Modules:
imgsdir <- file.path(figsdir,"SVG")
if (!dir.exists(imgsdir)) {
	# Create the directory.
	dir.create(imgsdir,recursive=TRUE)
} else {
	# Remove any existing figures.
	invisible({ file.remove(list.files(imgsdir,full.name=TRUE)) })
}

# all modules
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# Check that we are connected to Cytoscape.
cytoscapePing()

## Loop to create graphs:
message("\nCreating Cytoscape graphs!")
pbar <- txtProgressBar(max=length(module_list),style=3)
for (module in names(modules)){ 
	nodes <- modules[[module]]
	createCytoscapeGraph(netw_g, ppi_g, nodes, module, 
			     netwdir=netwdir,imgsdir=imgsdir)
	setTxtProgressBar(pbar, value = match(module,names(modules)))
}
close(pbar)

## When done, save Cytoscape session.
## NOTE: When on WSL, need to use Windows path format bc
## Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0(file_prefix,"Modules.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

message("\nDone!")
