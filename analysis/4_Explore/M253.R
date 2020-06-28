#!/usr/bin/env Rscript

#' ---
#' title: 
#' description:
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## 
#---------------------------------------------------------------------

## INPUT

## OPTION

## DEFAULT

## OUTPUT

#---------------------------------------------------------------------
## FUNCTIONS
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
## IMPORTS
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(geneLists)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

#---------------------------------------------------------------------
## LOAD DATA 
#---------------------------------------------------------------------

library(geneLists)
library(data.table)

data(iPSD)

data(gene_map)
data(partition)

# Map Uniprot to gene Symbols and Entrez IDs.
idx <- match(names(partition),gene_map$uniprot)
part_dt <- data.table(Module=partition,"Accession"=names(partition),
		      Symbol=gene_map$symbol[idx],
		      Entrez=gene_map$entrez[idx])
cb_module <- part_dt %>% filter(Symbol=="Arhgef9") %>% select(Module) %>% unlist()
part_dt$CbBioID <- part_dt$Entrez %in% iPSD$Arhgef9
part_dt %>% filter(Module == cb_module)

message("There are no Collybistin-BioID proteins in M253 (contains Arhgef9).")
message("Is this a false positive?")
