#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root = "~/projects/SwipProteomics"

## Options
n = "all" # the number of proteins to be analyzed or 'all'
save_rda = TRUE
drop_isSingle = TRUE # drop proteins with singular fits


## prepare the working environment ---------------------------------------------

# load renv
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})

## NOTE: my forks attempt to suppress much of MSstats verbosity as well as includes
## variants of the internal functions used by MSstatsTMT to fit protein-wise
## models and perform statistical comparisons between groups.


## load data ------------------------------------------------------------------

# load functions in root/R and make data in root/data accessible
devtools::load_all()

# load data in root/data
data(PD_raw) # 46 mb
data(gene_map)
data(PD_annotation)
data(msstats_contrasts)
# NOTE: msstats_contrasts is a matrix indicating pairwise contrasts between all
# BioFraction.Control and BioFraction.Mutant

## subset the data ------------------------------------------------------------
# NOTE: we have previously made sure Master.Protein.Accessions is a column of
# (single) Uniprot Identifiers. NOTE: Only run once!

# all proteins
proteins <- unique(as.character(PD_raw$Master.Protein.Accessions))

# define proteins to be analyzed:
if (n > length(proteins) | n == "all") {	
	# using all proteins
	message("\nAnalyzing 'all' (N=",
		formatC(length(proteins),big.mark=","),") proteins.")
	msstats_input <- PD_raw
} else {
	# subset the data
	message("\nAnalyzing a subset of data (n=",n," proteins).")
	prots <- sample(proteins,n)
	PD_filt <- PD_raw %>% 
		filter(Master.Protein.Accessions %in% prots)
	msstats_input <- PD_filt
} # EIS


## [1] convert to msstats format -----------------------------------------------
# MSstatsTMT key steps: 1-3.
# Proteins with a single feature are removed.

t0 = Sys.time()

suppressMessages({
PD_msstats <- PDtoMSstatsTMTFormat(msstats_input, 
				      PD_annotation, 
				      which.proteinid="Master.Protein.Accessions",
				      rmProtein_with1Feature = TRUE)
})

# aprox 7 minutes for all proteins
message("\nTime to pre-process ", n, " proteins: ", 
	round(difftime(Sys.time(), t0, units="min"),3)," minutes.")


## [2] summarize protein level data ----------------------------------------------

t0 = Sys.time()

suppressMessages({
	msstats_prot <- proteinSummarization(PD_msstats,
				     method="msstats",	
				     global_norm=TRUE,	
				     reference_norm=TRUE,
				     clusters=23)
})

# This takes about 11 minutes for 8.5 k proteins with 23 cores
message("\nTime to summarize ", n, " proteins: ", 
  round(difftime(Sys.time(),t0,units="min"),3)," minutes.")


## [3] perform statistical comparisons ----------------------------------------
# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# Abundance ~ (1|Mixture) + Condition

t0 = Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({
	msstats_results <- groupComparisonTMT(msstats_prot, msstats_contrasts)	
  })
})

# This takes about 21 minutes for 8.5 k proteins
message("\nTime to perform group comparisons for ", n, " proteins: ", 
	round(difftime(Sys.time(),t0,units="min"),3)," minutes.")


## annotate msstats_prot with gene ids ----------------------------------------

idx <- match(msstats_prot$Protein,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
msstats_prot <- tibble::add_column(msstats_prot, Symbol, .after="Protein")
msstats_prot <- tibble::add_column(msstats_prot, Entrez, .after="Symbol")


## combine normalized protein and statistical results -------------------------

if (drop_isSingle) {
	# if not NA, then issue. e.g. isSingleMeasure (7x in 100 protein fits)
	msstats_results <- msstats_results %>% filter(is.na(issue))
}

# drop issue column
msstats_results$issue <- NULL

# Annotated results and protein with BioFraction
msstats_results$BioFraction <- gsub("Mutant\\.F[0-9]{1,2}-Control\\.","",msstats_results$Label)
msstats_prot$BioFraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
msstats_prot$Genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)


## save results ---------------------------------------------------------------

if (save_rda) {

  myfile <- file.path(root,"data","msstats_results.rda")
  save(msstats_results,file=myfile,version=2)
  message("\nSaved ",basename(myfile),"in",dirname(myfile))
  
  # save normalized protein data in root/data
  myfile <- file.path(root,"data","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message("\nSaved ",basename(myfile),"in",dirname(myfile))
  
  # save the raw data in MSstatsTMT's format
  myfile <- file.path(root,"rdata","PD_msstats.rda")
  save(PD_msstats,file=myfile,version=2)
  message("\nSaved ",basename(myfile),"in",dirname(myfile))

}

message("\nCompleted MSstatsTMT intrafraction statistical analysis.")
