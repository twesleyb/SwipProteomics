#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Analysis options:

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

## THIS SCRIPT IS NOT COMPLETE. Manually processed network after 
## removing edges less than 0.9 and sending to cytoscape.

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
  library(RCy3) # For talking to Cytoscape.
  library(dplyr) # For manipulating data.
  library(igraph) # For creating graphs.
  library(data.table) # For working with tables.
  library(parallel) # For parallel proccessing.
  library(doParallel)
})

# Functions.
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

# Load the data from root/data.
data(tmt_protein)

# Load the graph partition:
data(partition)

# Load wash interactome.
data(wash_interactome)
wash_prots <- unique(wash_interactome$Accession) # Get uniprot accession

# Load networks.
data(ne_adjm) # loads "edges", then cast to adjm.
ne_adjm <- convert_to_adjm(edges)
data(ppi_adjm)
ppi_adjm <- convert_to_adjm(edges)

# Load gene map.
data(gene_map)

# Load module stats.
data(module_stats)

# Load sig prots.
data(sig_proteins)

#--------------------------------------------------------------------
## Create igraph graph objects.
#--------------------------------------------------------------------

# Create a list of all modules. 
module_list <- split(names(partition),partition)[-1] # drop M0
names(module_list) <- paste0("M",names(module_list))

# Insure that matrices are in matching order.
check <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!check) { stop() }

# Create igraph graph objects.
netw_g <- graph_from_adjacency_matrix(ne_adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)

ppi_g <- graph_from_adjacency_matrix(ppi_adjm,mode="undirected",diag=FALSE,
				     weighted=TRUE)

# Annotate graph's with gene symols.
symbols <- gene_map$"Gene Symbol"[match(names(V(netw_g)),gene_map$"Uniprot Accession")]
netw_g <- set_vertex_attr(netw_g,"symbol",value = symbols)
symbols <- gene_map$"Gene Symbol"[match(names(V(ppi_g)),gene_map$"Uniprot Accession")]
ppi_g <- set_vertex_attr(ppi_g,"symbol",value = symbols)

#--------------------------------------------------------------------
## Annotate graphs with additional meta data.
#--------------------------------------------------------------------

# Collect meta data from tmt_protein.
tmp_dt <- data.table(Accession = names(V(netw_g)),
		  Module = paste0("M",partition[names(V(netw_g))]))
noa <- left_join(tmp_dt, tmt_protein, by = "Accession") %>% 
	filter(!duplicated(Accession))
noa <- noa %>% select(Accession, Symbol, Entrez, Module, Adjusted.logFC, 
		      Adjusted.PercentWT, Adjusted.F, Adjusted.PValue, 
		      Adjusted.FDR)

# Add module colors.
data(module_colors)
noa$Color <- module_colors[noa$Module]

# Add WASH annotation.
noa$isWASH <- as.numeric(noa$Accession %in% wash_prots)

# Add sig prot annotations.
noa$sig85 <- as.numeric(noa$Accession %in% sig_proteins$sig85)
noa$sig62 <- as.numeric(noa$Accession %in% sig_proteins$sig62)
noa$sig968 <- as.numeric(noa$Accession %in% sig_proteins$sig968)

# Loop to add node attributes to netw_graph.
for (i in c(1:ncol(noa))) {
	namen <- colnames(noa)[i]
	col_data <- setNames(noa[[i]],nm=noa$Accession)
	netw_g <- set_vertex_attr(netw_g,namen,value=col_data[names(V(netw_g))])
}

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------

## Defaults.
n_cutoffs=10
netw_layout=NULL
netw_layout='force-directed edgeAttribute=weight' 

# Check that we are connected to Cytoscape.
response <- cytoscapePing()
if (response != "You are connected to Cytoscape!") { 
	stop("Start Cytoscape first!") 
}

# Prune weak edges.
# Seq from min(edge.weight) to max to generate cutoffs.
n_edges <- length(E(netw_g))
min_weight <- min(E(netw_g)$weight)
max_weight <- max(E(netw_g)$weight)
cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)

# Define a function that checks if graph is connnected at a 
# given edge weight threshold.
is_connected <- function(graph,threshold) {
	filt_graph <- delete.edges(graph, 
			       which(E(graph)$weight <= threshold))
	return(is.connected(filt_graph))
}

# Check if graph is connected or not at various thresholds.
# FIXME: Speed up by parallezing?
# NOTE: very slow for large graph.

## Why so slowwwww???
#nThreads <- parallel::detectCores()
#workers <- makeCluster(c(rep("localhost", nThreads)), type = "SOCK")
#registerDoParallel(workers)
#t0=Sys.time()
#checks <- foreach(i=cutoffs) %dopar% { is_connected(netw_g,i) }
#difftime(Sys.time(),t0) # Checking each threshold in serial takes about 6 seconds.
#stopCluster(workers)

# Limit is max(cutoff) at which the graph is still connected.
limit <- cutoffs[max(which(checks==TRUE))]
if (all(checks)) { stop("Error thesholding graph.") }

# Prune edges. NOTE: This removes all edge types.
limit = 0.9 # 500, 000 edges. is this possible.:
subg <- delete.edges(netw_g, which(E(netw_g)$weight <= limit))
n_edges <- length(E(subg))

# Write graph to file this is faster than sending to cytoscape.
myfile <- file.path(netwdir, paste0("network", ".gml"))
write_graph(subg, myfile, format = "gml")

# Send to Cytoscape.
# NOTE: underscores in attribute names are removed. 
winfile <- gsub("/mnt/d/", "D:/", myfile)
cysnetw <- importNetworkFromFile(winfile)
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
sig <- names(V(g))[V(g)$sig85 == 1 | V(g)$sig62 == 1 | V(g)$sig968 == 1]
ns <- names(V(g))[names(V(g)) %notin% sig]
if (length(ns) > 0) {
	setNodeColorBypass(ns,new.colors=col2hex("gray"))
}

# Bold border of BioID proteins.
#selectNodes(by="name",names(V(g))[V(g)$isWASH==1])
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

#--------------------------------------------------------------------------
## When done, save Cytoscape session.
#--------------------------------------------------------------------------

# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Modules.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

message("\nDone!")
