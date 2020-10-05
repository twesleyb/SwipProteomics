#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Options
save_rda = FALSE
n = 100

## ----------------------------------------------------------------------------
# prepare the working environment

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})


## ----------------------------------------------------------------------------
# load preprocessed data, annotation file, and contrasts

load(file.path(root,"rdata","msstats_input.rda"))
load(file.path(root,"rdata","msstats_annotation.rda"))
load(file.path(root,"rdata","msstats_contrasts.rda"))

# subset the data
# NOTE: we have previously made sure Master.Protein.Accessions is a column of
# (single) Uniprot Identifiers
# NOTE: Only run once!
proteins <- unique(as.character(msstats_input$Master.Protein.Accessions))
if (n > length(proteins) | n == "all") {	
	# using all proteins
	message("\nAnalyzing 'all' (N=",
		formatC(length(proteins),big.mark=","),") proteins.")
} else {
	# subset the data
	message("\nAnalyzing a subset of data (n=",n," proteins).")
	prots <- sample(proteins,n)
	msstats_input <- msstats_input %>% 
		filter(Master.Protein.Accessions %in% prots)
}

## ----------------------------------------------------------------------------
# convert to msstats format
t0 = Sys.time()
msstats_raw <- PDtoMSstatsTMTFormat(msstats_input, 
				    msstats_annotation, 
				    which.proteinid="Master.Protein.Accessions")	
message("Time to preprocess ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="s"),3)," seconds.")


## ----------------------------------------------------------------------------
# summarize protein level data
t0 = Sys.time()
suppressMessages({
msstats_prot <- proteinSummarization(msstats_raw,	
                                      method="msstats",	
                                      global_norm=TRUE,	
                                      reference_norm=TRUE)
})
message("Time to summarize ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="s"),3)," seconds.")

## ----------------------------------------------------------------------------
# perform statistical comparisons

t0 = Sys.time()
suppressWarnings({
	suppressMessages({
  msstats_results <- groupComparisonTMT(msstats_prot, msstats_contrasts)	
	})
})
message("Time to perform group comparisons for ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="s"),3)," seconds.")


## ----------------------------------------------------------------------------
# save results

if (save_rda) {

  myfile <- file.path(root,"rdata","msstats_results.rda")
  save(msstats_results,file=myfile,version=2)
  message("Saved",basename(myfile),"in",dirname(myfile))
  
  myfile <- file.path(root,"rdata","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message("Saved",basename(myfile),"in",dirname(myfile))
  
  myfile <- file.path(root,"rdata","msstats_raw.rda")
  save(msstats_raw,file=myfile,version=2)
  message("Saved",basename(myfile),"in",dirname(myfile))

}

message("\nCompleted MSstatsTMT intrafraction statistical analysis.")
