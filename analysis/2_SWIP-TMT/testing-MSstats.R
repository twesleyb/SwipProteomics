#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# data in root/rdata:
# * data_prot.rda

## Prepare the working environment --------------------------------------------

root <- "~/projects/SwipProteomics"
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata") # temporary/large files not tracked by git


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(MSstatsTMT)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## Perform protein-level statistical tests --------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment.

# load the preprocessed MSstats protein level data
myfile <- file.path(datadir,"msstats_prot.rda")
load(myfile); data_prot <- msstats_prot

# test for all the possible pairs of conditions	
#all_results <- groupComparisonTMT(data_prot)	
# FIXME: need to pass contrast matrix!

myfile <- file.path(datadir,"msstats_prot.rda")
load(myfile); data_prot <- msstats_prot
data = data_prot
contrast.matrix = 'pairwise'
moderated = FALSE
adj.method = 'BH'
remove_norm_channel = TRUE
remove_empty_channel = TRUE

###############################################################################
# Work through statistics:
# OBJECTIVE: implement mixed effect model for repeated measures:
# fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=data2) # [2]
