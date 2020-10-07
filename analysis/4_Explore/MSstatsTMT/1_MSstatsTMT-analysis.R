#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Options
n = "all"
save_rda = TRUE

## prepare the working environment ---------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})


## load preprocessed data, annotation file, and contrasts ---------------------

load(file.path(root,"rdata","PD_input.rda"))
load(file.path(root,"rdata","PD_annotation.rda"))
load(file.path(root,"rdata","msstats_contrasts.rda"))

## subset the data ------------------------------------------------------------
# NOTE: we have previously made sure Master.Protein.Accessions is a column of
# (single) Uniprot Identifiers
# NOTE: Only run once!

# all proteins
proteins <- unique(as.character(PD_input$Master.Protein.Accessions))

# define proteins to be analyzed:
if (n > length(proteins) | n == "all") {	
	# using all proteins
	message("\nAnalyzing 'all' (N=",
		formatC(length(proteins),big.mark=","),") proteins.")
} else {
	# subset the data
	message("\nAnalyzing a subset of data (n=",n," proteins).")
	prots <- sample(proteins,n)
	PD_input <- PD_input %>% 
		filter(Master.Protein.Accessions %in% prots)
} # EIS


## [1] convert to msstats format -----------------------------------------------
# MSstatsTMT key steps: 1-3.

# This step takes about 7 minutes for 8.5 k proteins
t0 = Sys.time()

msstats_input <- PDtoMSstatsTMTFormat(PD_input, 
			    PD_annotation, 
			    which.proteinid="Master.Protein.Accessions")	

message("Time to preprocess ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="s"),3)," seconds.")


## [2] summarize protein level data ----------------------------------------------

# This takes about 45 minutes for 8.5 k proteins
t0 = Sys.time()

suppressMessages({
	msstats_prot <- proteinSummarization(msstats_input,
					     method="msstats",	
					     global_norm=TRUE,	
					     reference_norm=TRUE)
})

message("Time to summarize ", n, " proteins: ", 
  round(difftime(Sys.time(),t0,units="s"),3)," seconds.")


## [3] perform statistical comparisons ----------------------------------------
# This takes about 21 minutes for 8.5 k proteins

t0 = Sys.time()

suppressWarnings({
	suppressMessages({
  msstats_results <- groupComparisonTMT(msstats_prot, msstats_contrasts)	
	})
})

message("Time to perform group comparisons for ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="s"),3)," seconds.")


## save results ---------------------------------------------------------------

if (save_rda) {
  # save output as rda in root/rdata
  myfile <- file.path(root,"rdata","msstats_results.rda")
  save(msstats_results,file=myfile,version=2)
  message("Saved ",basename(myfile),"in",dirname(myfile))
  
  myfile <- file.path(root,"rdata","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message("Saved ",basename(myfile),"in",dirname(myfile))
  
  myfile <- file.path(root,"rdata","msstats_input.rda")
  save(msstats_input,file=myfile,version=2)
  message("Saved ",basename(myfile),"in",dirname(myfile))

}

message("\nCompleted MSstatsTMT intrafraction statistical analysis.")
