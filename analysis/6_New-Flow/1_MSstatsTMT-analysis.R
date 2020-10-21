#!/usr/bin/env Rscript

# title: MSstatsTMT 
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root = "~/projects/SwipProteomics"

## Options
nprot = "all" # the number of proteins to be analyzed or 'all'
FDR_alpha = 0.05 # FDR threshold for significance
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


## NOTE: my fork attempts to suppress much of MSstats's verbosity as well as 
## includes access to the internal functions used by MSstatsTMT to fit 
## protein-wise models and perform statistical comparisons between groups.
## MSstatsTMT is a wrapper around MSstats. My fork allows you to pass
## arguments for parallel processing to proteinSummarization to speed things up.


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
	message("\nAnalyzing a subset of data (n=",nprot,") proteins.")
	prots <- sample(proteins,nprot)
	pd_filt <- pd_psm %>% 
		filter(Master.Protein.Accessions %in% prots)
	msstats_input <- pd_filt
} # EIS


## [1] convert to msstats format -----------------------------------------------
# Proteins with a single feature are removed.

t0 = Sys.time()

suppressMessages({ # verbosity

  msstats_psm <- PDtoMSstatsTMTFormat(msstats_input, 
				      pd_annotation, 
				      which.proteinid="Master.Protein.Accessions",
				      rmProtein_with1Feature = TRUE)
})

# aprox 7 minutes for all proteins
message("\nTime to pre-process ", nprot, " proteins: ", 
	round(difftime(Sys.time(), t0, units="min"),3)," minutes.")


## [2] summarize protein level data ----------------------------------------------
# Perform protein summarization for each run.
t0 = Sys.time()

suppressMessages({ # verbosity

  msstats_prot <- proteinSummarization(msstats_psm,
				       method="msstats",	
				       global_norm=TRUE,	
				       MBimpute=TRUE,
				       reference_norm=TRUE,
				       clusters=23)

})

# This takes about 11 minutes for 8.5 k proteins with 23 cores
# FIXME: fix warnings messages about closing clusters.
message("\nTime to summarize ", nprot, " proteins: ", 
  round(difftime(Sys.time(), t0, units="min"), 3)," minutes.")


## [3] perform statistical comparisons ----------------------------------------
# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# lmerTest::lmer(Abundance ~ (1|Mixture) + Condition)
# We specify Condition as Genotype.BioFraction for all intra-fraction
# comparisons. T-statistics are moderated using ebayes methods in limma.

t0 = Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity

    msstats_results <- groupComparisonTMT(msstats_prot, 
					  msstats_contrasts, 
					  moderated = TRUE)

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


## clean-up msstats_results ---------------------------------------------------

# if not NA, then issue. e.g. isSingleMeasure (7x in 100 protein fits)
msstats_results <- msstats_results %>% filter(is.na(issue))

# drop issue column
msstats_results$issue <- NULL

# annotate results with gene symbols
idx <- match(msstats_results$Protein,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
msstats_results <- tibble::add_column(msstats_results,Symbol,.after="Protein")

# summary
message("\nSummary of signifcant (FDR<",
	FDR_alpha,") proteins for intrafraction comparisons:")
df = msstats_results %>% 
	group_by(Label) %>% arrange(adj.pvalue) %>% 
	summarize(`n Sig` = sum(adj.pvalue < FDR_alpha),
		  `Top 5 Sig Prots` = paste(head(Symbol[adj.pvalue < FDR_alpha]),
					    collapse=", "), .groups="drop")
colnames(df)[1] <- "Contrast"
knitr::kable(df)


## save results ---------------------------------------------------------------

if (save_rda) {

  # save msstats_results -- MSstatsTMT statistical results
  myfile <- file.path(root,"data","msstats_results.rda")
  save(msstats_results,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))
  
  # save msstats_prot -- the normalized protein
  myfile <- file.path(root,"data","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))
  
  # save msstats_psm -- the psm level data reformatted for MSstats
  myfile <- file.path(root,"rdata","msstats_psm.rda")
  save(msstats_psm,file=myfile,version=2)
  message("\nSaved ",basename(myfile)," in ",dirname(myfile))

}

message("\nCompleted MSstatsTMT intrafraction statistical analysis.")
