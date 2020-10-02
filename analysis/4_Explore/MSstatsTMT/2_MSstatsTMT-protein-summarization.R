#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

# FIXME: catch/suppress warnings about singular fit?
# FIXME: exchange pbar for messages

## Input ----------------------------------------------------------------------
# specify the projects root directory:
root = "~/projects/SwipProteomics"

# MSstatsTMT preprocessed data:
# * data_pd.rda in root/rdata

## Options --------------------------------------------------------------------
save_rda = TRUE # save key R objects?

## Output ---------------------------------------------------------------------
# * several intermediate datasets in root/rdata
# * the MSstatsTMT statistical results for intrafraction comparisons saved as
#   an excel workbook in root/tables

## FUNCTIONS -----------------------------------------------------------------

squote <- function(string) { 
	# wrap a string in singular quotation marks (')
	return(paste0("'",string,"'"))
}


mkdir <- function(...,warn=TRUE,report=FALSE) {
	# create a new directory
	newdir <- file.path(...)
	if (warn & dir.exists(newdir)) {
		warning("dir exists")
	} else if (!dir.exists(newdir)) { 
		dir.create(newdir)
		if (report) {
		message(paste("created",newdir))
		}
	}
}


## Prepare the working environment --------------------------------------------

# project directories
datadir <- file.path(root,"data"); mkdir(datadir,warn=FALSE)
rdatdir <- file.path(root,"rdata"); mkdir(rdatdir,warn=FALSE)
downdir <- file.path(root,"downloads"); mkdir(downdir,warn=FALSE)


# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	suppressWarnings({ library(getPPIs) }) # FIXME: annoying warnings!
	library(data.table)
	library(MSstatsTMT) # twesleyb/MSstatsTMT
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## load preprocessed data ------------------------------------------------------

load(file.path(root,"rdata","data_pd.rda"))

# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
msstats_prot <- proteinSummarization(data_pd,
			  method="msstats",	
			  global_norm=TRUE,	
			  reference_norm=TRUE,	
			  remove_norm_channel = TRUE)

if (save_rda) {
  # save to file
  myfile <- file.path(rdatdir, "msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message(paste("\nSaved",myfile))
}
