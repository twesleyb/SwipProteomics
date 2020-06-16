#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
  library(RCy3)
  library(dplyr)
  library(igraph)
  library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

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
data(adjm)
data(ne_adjm)
data(ppi_adjm)

# Load gene map
gene_map <- readRDS(file.path(rdatdir,"gene_map.RData"))

#--------------------------------------------------------------------
## Create igraph graph objects.
#--------------------------------------------------------------------

# List of all modules. 
module_list <- split(names(partition),partition)[-1] # drop M0
names(module_list) <- paste0("M",names(module_list))

# Insure that matrices are in matching order.
col_names <- colnames(adjm)
ne_adjm <- ne_adjm[col_names,col_names]
ppi_adjm <- ppi_adjm[col_names,col_names]

# Create igraph graphs from adjacency matrices.
adjm_g <- graph_from_adjacency_matrix(adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
netw_g <- graph_from_adjacency_matrix(ne_adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
ppi_g <- graph_from_adjacency_matrix(ppi_adjm,mode="undirected",diag=FALSE,
				     weighted=TRUE)

# Annotate graph's with gene symols.
symbols <- gene_map$symbol[match(names(V(adjm_g)),gene_map$uniprot)]
adjm_g <- set_vertex_attr(adjm_g,"symbol",value = symbols)
symbols <- gene_map$symbol[match(names(V(netw_g)),gene_map$uniprot)]
netw_g <- set_vertex_attr(netw_g,"symbol",value = symbols)
symbols <- gene_map$symbol[match(names(V(ppi_g)),gene_map$uniprot)]
ppi_g <- set_vertex_attr(ppi_g,"symbol",value = symbols)

#--------------------------------------------------------------------
## Annotate graphs with additional meta data.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------

## Loop to create graphs:
netwdir = file.path(root,"networks")
for (module_name in names(module_list)){
	nodes = module_list[[module_name]]
	createCytoscapeGraph(netw_g,ppi_g,nodes,module_name,netwdir=netwdir)
}

# When done, save cytoscape session.
myfile <- file.path(netwdir,paste0("Modules.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile)
saveSession(winfile)
