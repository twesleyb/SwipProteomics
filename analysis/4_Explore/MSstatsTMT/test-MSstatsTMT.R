#!/usr/bin/env Rscript

# title: MSstatsTMT 
# NOTE: adapted from MSstatsTMT vignette 

root <- "~/projects/SwipProteomics"
renv::load(root)


library(MSstatsTMT)	

data(raw.pd)
data(annotation.pd)
	
# 1.
input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)	

# 2.
quant.msstats <- proteinSummarization(input.pd,	
                                      method="msstats",	
                                      global_norm=TRUE,	
                                      reference_norm=TRUE,	
                                      remove_norm_channel = TRUE,	
                                      remove_empty_channel = TRUE)	
	
# 3.
test.pairwise <- groupComparisonTMT(quant.msstats)	
