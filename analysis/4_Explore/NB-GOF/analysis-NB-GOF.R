#!/usr/bin/env Rscript

# title: SwipProteomics
# description: anaysis of the goodness of fit of negative binomial models
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
#root <- Sys.getEnv("ROOT") # set by activate script
root <- getrd()
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(NBGOF)
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# load data from NBGOF package
data(arab)
data(earthquake)
