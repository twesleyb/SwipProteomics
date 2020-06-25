#!/usr/bin/env bash

#--------------------------------------------------------------------
## Misc function - getrd
#--------------------------------------------------------------------
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

#--------------------------------------------------------------------
## Preepare the workspace.
#--------------------------------------------------------------------

root <- getrd()
renv::load(root,quiet=TRUE)

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

devtools::load_all()

# load the data
data(tmt_protein)

# Plot volcano plot.
plot <- ggplot()
plot <- ggplot(tmt_protein, aes(x = Adjusted.logFC, y = -log10(Adjusted.PValue)))
plot <- plot + geom_point()





