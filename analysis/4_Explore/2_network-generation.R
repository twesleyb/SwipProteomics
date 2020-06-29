#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: generate protein networks
#' authors: Tyler W A Bradshaw
#' ---

## INPUT:
# tmt_protein in root/data

## OPTIONS:

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# Get the repository's root directory.
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
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(WGCNA)
	library(neten)
	library(getPPIs)
	library(data.table)
})

# Project imports.
devtools::load_all()

#--------------------------------------------------------------------
## Create protein covariation network.
#--------------------------------------------------------------------

# Load the proteomics data.
data(tmt_protein)

# Cast to a data.matrix.
dm <- tmt_protein %>% as.data.table() %>%
	dcast(Sample ~ Accession, value.var="Intensity") %>%
	as.matrix(rownames=TRUE) %>% log2()

# Create correlation (adjacency) matrix.
message("\nGenerating protein co-variation matrix using bicor().")
adjm <- WGCNA::bicor(dm)

# Enhanced network.
message("\nPerforming network enhancement with to denoise network.")
ne_adjm <- neten::neten(adjm)

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

message("\nSaving the data.")

# Save adjm as csv.
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","adjm.csv"))

# Save enhanced adjm as csv.
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ne_adjm.csv"))
