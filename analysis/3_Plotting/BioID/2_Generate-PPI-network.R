#!/usr/bin/env Rscript

# title: WASH iBioID Proteomics Analysis
# description: building a PPI network for WASH proteome.
# author: Tyler W A Bradshaw

## User parameters to change:
FDR_alpha = 0.1
enrichment_threshold = log2(3.0)

## Input data in root/rdata
# * bioid_se

## Output in root/manuscript/tables:
# * WASH_Network_PPIs.xlsx

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(igraph) # For networks.
	library(getPPIs) # For mapping gene names.
	library(data.table) # For data.tables.
	library(openxlsx) # For working with .xlsx files.
})

# Directories.
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"manuscript","tables")

# Load any additional functions in root/R.
devtools::load_all()

# Load BioID results.
data(bioid_se)
data(bioid_gene_map) # gene_map
results <- DEP::get_results(bioid_se$results)


#---------------------------------------------------------------------
## Create PPI graph.
#---------------------------------------------------------------------

# Get entrez identifiers for proteins of interest.
sig <- results$WASH_vs_Control_p.adj < FDR_alpha
up <- results$WASH_vs_Control_ratio > enrichment_threshold
uniprot <- results$ID[sig & up]
entrez <- gene_map$entrez[match(uniprot,gene_map$uniprot)]

# Load PPIs from getPPIs.
data(musInteractome)

# Subset dataset--get PPIs among our proteins of entrest.
ppis <- musInteractome %>% 
	filter(osEntrezA %in% entrez & osEntrezB %in% entrez) %>% 
	select(osEntrezA,osEntrezB,Interactor_A_Taxonomy,
	       Interaction_detection_methods,
	       Publications,Interaction_type,Confidence_score,
	       Source_database,Methods)

# Create edges data.frame -- aka SIF or simple interaction file.
sif <- ppis 
colnames(sif)[c(1,2,3)] <- c("EntrezA","EntrezB","Taxonomy")
sif$Weight <- 1

# Map Entrez to protein ids.
sif$ProteinA <- gene_map$symbol[match(sif$EntrezA,gene_map$entrez)]
sif$ProteinB <- gene_map$symbol[match(sif$EntrezB,gene_map$entrez)]

# Node attributes -- first column must be matching protein identifiers.
noa <- data.table(entrez = entrez, protein = uniprot, 
		  symbol=results$name[sig & up])
rownames(noa) <- noa$entrez

# Create graph.
message("\nCreating PPI graph from WASH iBioID interactome.")
g <- graph_from_data_frame(sif,vertices=noa,directed=FALSE)

# Summarize number of nodes and edges.
message(paste("\nNumber of Nodes:",length(V(g))))
message(paste("\nNumber of Edges:",length(E(g))))
