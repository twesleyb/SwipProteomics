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

# Load NDD associated proteins.
data(NDD_proteins)

# Load networks.
data(ne_adjm) # loads "edges", then cast to adjm.
ne_adjm <- convert_to_adjm(edges)
data(ppi_adjm)
ppi_adjm <- convert_to_adjm(edges)

# Load gene map.
data(gene_map)

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
symbols <- gene_map$symbol[match(names(V(netw_g)),gene_map$uniprot)]
netw_g <- set_vertex_attr(netw_g,"symbol",value = symbols)
symbols <- gene_map$symbol[match(names(V(ppi_g)),gene_map$uniprot)]
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

# Add NDD annotations.
noa$isNDD <- as.numeric(noa$Accession %in% names(NDD_proteins))
noa$NDD <- NDD_proteins[noa$Accession]
noa$NDD[is.na(noa$NDD)] <- "none"

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

# Network images will be saved in networks/Modules:
imgsdir <- file.path(figsdir,"SVG")
if (!dir.exists(imgsdir)) {
	# Create the directory.
	dir.create(imgsdir)
} else {
	# Remove any existing figures.
	invisible({ file.remove(list.files(imgsdir,full.name=TRUE)) })
}

# Loop to create graphs:
message("\nCreating Cytoscape graphs!")
pbar <- txtProgressBar(max=length(module_list),style=3)
for (module_name in names(module_list)){ 
	nodes <- module_list[[module_name]]
	createCytoscapeGraph(netw_g, ppi_g, nodes, module_name, 
			     netwdir=netwdir,imgsdir=imgsdir)
	setTxtProgressBar(pbar, value = match(module_name,names(module_list)))
}
close(pbar)

# When done, save Cytoscape session.
# NOTE: When on WSL, need to use Windows path format bc
# Cytoscape is a Windows program.
myfile <- file.path(netwdir,paste0("Modules.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile) 
saveSession(winfile)

message("\nDone!")
