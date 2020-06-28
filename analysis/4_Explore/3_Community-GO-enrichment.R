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

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

# Load the partition and tmt data.
data(partition)
data(tmt_protein)
data(sig_modules)
data(module_stats)

myfile <- file.path(rdatdir,"SigModule_partition.csv")
module_partition <- fread(myfile,drop=1) %>% unlist() + 1

# Annotate modules as either up or down.
data(module_stats)
module_stats$Module <- paste0("M",module_stats$Module)
substats <- module_stats %>% filter(module_stats$PAdjust < 0.05)
substats$down <- substats$logFC < 0
substats$Module_Group <- module_partition[substats$Module]
substats <- substats %>% arrange(Module_Group)

modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# collect list of entrez ids and gene symbols.
gene_list <- lapply(modules,function(x) {
			    tmt_protein$Entrez[match(x,tmt_protein$Accession)]
			    })

#  Split data by module group and up/down.
data_list <- substats %>% group_by(Module_Group,down) %>% group_split()

gene_list <- sapply(data_list,function(x) unlist(gene_list[x$Module]))

updown = c("TRUE"="down","FALSE"="up")
names(gene_list) <- sapply(data_list,function(x) paste0("G",unique(x$Module_Group),":",updown[as.character(unique(x$down))]))

# Build a GO reference collection:
message("\nBuilding a mouse GO reference collection with anRichment.")
gocollection <- suppressPackageStartupMessages({
	anRichment::buildGOcollection(organism="mouse")
})

# perform module go enrichment.
# NOTE: background defaults to all genes in gene_list.
go_gse <- gse(gene_list,gocollection)

# Save data.
mydata <- bind_rows(go_gse,.id="Group")
fwrite(mydata,file.path(rdatdir,"community_GO_results.csv")
