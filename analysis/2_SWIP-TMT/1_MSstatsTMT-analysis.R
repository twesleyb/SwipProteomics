#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root = "~/projects/SwipProteomics"

## Options
<<<<<<< HEAD
nprot = "all" # the number of proteins to be analyzed or 'all'
=======
nprot = 100 #"all" # the number of proteins to be analyzed or 'all'
>>>>>>> dev
save_rda = TRUE


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
suppressWarnings({ devtools::load_all() })

# load data in root/data
data(pd_psm) # 46 mb
data(gene_map)
data(pd_annotation)
data(msstats_contrasts)
# NOTE: msstats_contrasts is a matrix indicating pairwise contrasts between all
# BioFraction.Control and BioFraction.Mutant

## subset the data ------------------------------------------------------------
# NOTE: we have previously made sure Master.Protein.Accessions is a column of
# (single) Uniprot Identifiers. NOTE: Only run once!

# all proteins
proteins <- unique(as.character(pd_psm$Master.Protein.Accessions))

# define proteins to be analyzed:
if (nprot > length(proteins) | nprot == "all") {	
	# using all proteins
	message("\nAnalyzing 'all' (N=",
		formatC(length(proteins),big.mark=","),") proteins.")
	msstats_input <- pd_psm
} else {
	# subset the data
<<<<<<< HEAD
	message("\nAnalyzing a subset of data (n=",nprot," proteins).")
	prots <- sample(proteins,nprot)
	PD_filt <- pd_psm %>% 
=======
	message("\nAnalyzing a subset of data (n=",nprot,") proteins.")
	prots <- sample(proteins,nprot)
	pd_filt <- pd_psm %>% 
>>>>>>> dev
		filter(Master.Protein.Accessions %in% prots)
	msstats_input <- pd_filt
} # EIS


## [1] convert to msstats format -----------------------------------------------
# MSstatsTMT key steps: 1-3.
# Proteins with a single feature are removed.

t0 = Sys.time()

suppressMessages({
<<<<<<< HEAD
msstats_psm <- PDtoMSstatsTMTFormat(pd_psm, 
				    pd_annotation, 
				    which.proteinid="Master.Protein.Accessions",
				    rmProtein_with1Feature = TRUE)
=======
msstats_psm <- PDtoMSstatsTMTFormat(msstats_input, 
				   pd_annotation, 
				   which.proteinid="Master.Protein.Accessions",
				   rmProtein_with1Feature = TRUE)
>>>>>>> dev
})

# aprox 7 minutes for all proteins
message("\nTime to pre-process ", nprot, " proteins: ", 
	round(difftime(Sys.time(), t0, units="min"),3)," minutes.")


## [2] summarize protein level data ----------------------------------------------
# Perform protein summarization for each run.
t0 = Sys.time()

suppressMessages({
	msstats_prot <- proteinSummarization(msstats_psm,
				     method="msstats",	
				     global_norm=TRUE,	
				     reference_norm=TRUE,
				     clusters=23)
})

# This takes about 11 minutes for 8.5 k proteins with 23 cores
#FIXME: fix warnings messages about closing clusters
message("\nTime to summarize ", nprot, " proteins: ", 
  round(difftime(Sys.time(), t0, units="min"), 3)," minutes.")


## [3] perform statistical comparisons ----------------------------------------
# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# lmerTest::lmer(formula(Abundance ~ (1|Mixture) + Condition))

t0 = Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({
	msstats_results <- groupComparisonTMT(msstats_prot, msstats_contrasts)	
  })
})

# This takes about 21 minutes for 8.5 k proteins
message("\nTime to perform group comparisons for ", nprot, " proteins: ", 
	round(difftime(Sys.time(),t0 ,units="min"),3)," minutes.")


## annotate msstats_prot with gene ids ----------------------------------------

idx <- match(msstats_prot$Protein, gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
msstats_prot <- tibble::add_column(msstats_prot, Symbol, .after="Protein")
msstats_prot <- tibble::add_column(msstats_prot, Entrez, .after="Symbol")


## combine normalized protein and statistical results -------------------------

# if not NA, then issue. e.g. isSingleMeasure (7x in 100 protein fits)
msstats_results <- msstats_results %>% filter(is.na(issue))

# drop issue column
msstats_results$issue <- NULL

# Annotated results and protein with BioFraction
msstats_results$BioFraction <- gsub("Mutant\\.F[0-9]{1,2}-Control\\.","",
				    msstats_results$Label)
condition <- as.character(msstats_prot$Condition)
msstats_prot$BioFraction <- sapply(strsplit(condition,"\\."),"[",2)
msstats_prot$Genotype <- sapply(strsplit(condition,"\\."),"[",1)

<<<<<<< HEAD

## Calculate Protein abundance adjusted for fraction differences --------------

# NOTE: these adjusted values are not used for modeling, only 
# plotting purposes.

## WORK

# Calculate protein abundance, adjusted for fraction differences.
# FIXME: Is CPM step needed? Log is necessary, but CPM norm?
#logCPM <- edgeR::cpm(glm_results$dge, log=TRUE)

# Remove effect of fraction.
#dm <- limma::removeBatchEffect(logCPM,
#			       batch=dge$samples$fraction,
#			       design=model.matrix(~treatment,data=dge$samples))

# Collect results.
#dt <- as.data.table(dm,keep.rownames="Accession") %>% 
#		reshape2::melt(id.var="Accession",
#			       variable.name="Sample",
#			       value.name="Adjusted.Intensity")
#dt$Sample <- as.character(dt$Sample)

# Combine Adjusted protein with additional sample meta data.
#tmp_prot <- tmt_protein %>% filter(Treatment != "SPQC") %>% 
#	as.data.table()
#adjusted_prot <- left_join(tmp_prot,dt,
#			   by=intersect(colnames(tmp_prot),colnames(dt)))

# Add stats to adjusted protein data.
#adjusted_prot <- left_join(adjusted_prot,alt_results,
#			   by=intersect(colnames(adjusted_prot),
#					colnames(alt_results)))

# Unlog the data.
#adjusted_prot$Adjusted.Intensity <- 2^adjusted_prot$Adjusted.Intensity

=======
>>>>>>> dev

## save results ---------------------------------------------------------------

if (save_rda) {

  myfile <- file.path(root,"data","msstats_results.rda")
  save(msstats_results,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))
  
  # save normalized protein data in root/data
  myfile <- file.path(root,"data","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))
  
  # save the raw data in MSstatsTMT's format
<<<<<<< HEAD
  # NOTE: its too large to store in data and track with git.
  # save in root/rdata
  myfile <- file.path(root,"rdata","msstats_psm.rda")
  save(msstats_psm,file=myfile,version=2)
  message("\nSaved ",basename(myfile),"in",dirname(myfile))
=======
  myfile <- file.path(root,"rdata","msstats_psm.rda")
  save(msstats_psm,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))
>>>>>>> dev

}

message("\nCompleted MSstatsTMT intrafraction statistical analysis.")
