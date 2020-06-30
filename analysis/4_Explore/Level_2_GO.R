#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT proteomics
#' description: Analysis of modules for GO enrichment
#' authors: Tyler W A Bradshaw
#' ---

## User parameters to change:
FDR_alpha = 0.05 # significance threshold for gse enrichment.

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

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Load additional functions in root/R.
suppressWarnings({ devtools::load_all() })

# Directories.
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

# Load the partition and tmt data.
data(partition)
data(tmt_protein)
data(sig_proteins)
data(wash_interactome)

# Get list of up and down regulated proteins.
sig_prots <- sig_proteins
tmt_protein$isUp <- tmt_protein$Adjusted.logFC > 0 
data_list <- tmt_protein %>% 
	filter(Accession %in% sig_proteins) %>% 
	group_by(isUp) %>% group_split()
gene_list <- lapply(data_list,function(x) unique(x$Entrez))
names(gene_list) <- c("SigDown","SigUp")

#---------------------------------------------------------------------
# how many.
#---------------------------------------------------------------------

wash_prots <- wash_interactome$Accession

sig85 <- tmt_protein %>% filter(FDR<0.1) %>% select(Accession) %>% 
	unlist() %>% unique()
sum(wash_prots %in% sig85)
sum(wash_prots %in% sig_proteins)

#---------------------------------------------------------------------
## Module GO enrichment.
#---------------------------------------------------------------------

# Build a GO reference collection:
message("\nBuilding a mouse GO reference collection with anRichment.")
gocollection <- suppressPackageStartupMessages({
	anRichment::buildGOcollection(organism="mouse")
})

# perform module go enrichment.
# NOTE: background defaults to all genes in gene_list.
go_gse <- gse(gene_list,gocollection)

# save significant results.
write_excel(go_gse,file.path(tabsdir,"Level2_Up_Down_GO_Results.xlsx"))
