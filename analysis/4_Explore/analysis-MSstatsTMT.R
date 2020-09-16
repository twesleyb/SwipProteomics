#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT with MSstats
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Misc functions -------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# get the repository's root directory
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

# Prepare the R environment ---------------------------------------------------

# load renv
#root <- Sys.getEnv("ROOT") # set by activate script in 
root <- getrd()
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(MSstats)
	library(MSstatsTMT)
})

# load functions in root/R and data in root/data
#FIXME: where is 'Error in readRDS(file) : unknown input format' coming from?
suppressWarnings({ devtools::load_all() })

data(tidy_peptide) # tmt_protein
data(input.pd) # MSstatsTMT dataset

knitr::kable(str(input.pd))
knitr::kable(str(tidy_peptide))
