#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})

# load preprocessed data, annotation file, and contrasts
load(file.path(root,"rdata","msstats_input.rda"))
load(file.path(root,"rdata","mstats_annotation.rda"))
load(file.path(root,"rdata","contrast_matrix.rda"))
	
# convert to msstats format
msstats_raw <- PDtoMSstatsTMTFormat(msstats_input, msstats_annotation)	

# summarize protein level data
msstast_prot <- proteinSummarization(msstast_raw,	
                                      method="msstats",	
                                      global_norm=TRUE,	
                                      reference_norm=TRUE)
	
# perform statistical comparisons
msstats_results <- groupComparisonTMT(msstats_prot,contrast_matrix)	

# save
myfile <- file.path(root,"rdata","msstats_results.rda")
save(msstats_results,file=myfile,version=2)

myfile <- file.path(root,"rdata","msstats_prot.rda")
save(msstats_prot,file=myfile,version=2)

myfile <- file.path(root,"rdata","msstats_raw.rda")
save(msstats_raw,file=myfile,version=2)
