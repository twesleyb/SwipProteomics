#!/usr/bin/env Rscript

# title: MSstatsTMT 

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
library(MSstats) # twesleyb/MSstats
library(MSstatsTMT) # twesleyb/MSstats

# load preprocessed data, annotation file, and contrasts
load(file.path(root,"rdata","msstats_input.rda"))
load(file.path(root,"rdata","mstats_annotation.rda"))
load(file.path(root,"rdata","contrast_matrix.rda"))
	
msstats_raw <- PDtoMSstatsTMTFormat(msstats_input, msstats_annotation)	

msstast_prot <- proteinSummarization(msstast_raw,	
                                      method="msstats",	
                                      global_norm=TRUE,	
                                      reference_norm=TRUE)
	
msstats_resulst <- groupComparisonTMT(msstats_prot,contrast_matrix)	
