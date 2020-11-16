#!/usr/bin/env Rscript 

# title: MSstatsTMT analysis of SWIP-TMT Proteomics
# description: analysis of protein-level comparisons with MSstatsTMT
# author: twab

## MSstatsTMT options
MBimpute <- TRUE
rm_single <- TRUE
global_norm <- TRUE
reference_norm <- TRUE
remove_norm_channel <- TRUE

## threshold for significance
FDR_alpha = 0.05


## prepare the working environment ---------------------------------------------

# load functions in root/R and make data in root/data accessible
library(SwipProteomics)

# load data in root/data
data(pd_psm) 
data(pd_annotation)
data(mut_vs_control) # 'Mutant-Control' comparison
data(msstats_contrasts) # 'intra-BioFraction' comparisons

# NOTE: msstats_contrasts is a matrix specifying pairwise contrasts between all
# 'BioFraction.Control' and 'BioFraction.Mutant' Conditions.

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  # my MSstats forks:
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})

## NOTE: my fork attempts to remove much of MSstats's verbosity.
## You may find some debug.log files on your computer still, their source
## has not been found yet.
## twesleyb/MSstatsTMT includes access to the internal functions used by
## MSstatsTMT to fit protein-wise models and perform statistical comparisons.
## MSstatsTMT is a wrapper around MSstats. My fork allows you to pass arguments
## for parallel processing to MSstats::proteinSummarization to speed things up.


## [1] convert to msstats format -----------------------------------------------
# Proteins with a single feature (i.e. peptide) are removed if rm_single.

message("\nConverting PD PSM-level data into MSstatsTMT format.")

t0 <- Sys.time()

# aprox 7 minutes for all proteins
suppressMessages({ # verbosity
  msstats_psm <- PDtoMSstatsTMTFormat(pd_psm,
    pd_annotation, # see 0_PD-data-preprocess.R
    which.proteinid = "Master.Protein.Accessions",
    rmProtein_with1Feature = rm_single
  )
})

message("\nTime to reformat PSM data: ", 
	round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


## [added] QC-based PSM filtering ------------------------------------------------

# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated.
# These ratios are then binned based on mean(log2(Intensity)) into
# 5 bins. For each bin, measuremnts that are outside
# +/- nSD x standard deviations from the bin's mean are removed as outliers.

# This is essential because our normalization strategy depends upon the SPQC
# samples. We cannot perform normalization to SPQC PSM that are highly variable.

# identify psm outliers for each Mixture
mix <- c("M1","M2","M3")
outlier_list <- lapply(mix, rmOutlierPSM, msstats_psm, nbins=5, nSD=4)
names(outlier_list) <- mix

# all psm outliers
psm_outliers <- sapply(outlier_list,function(x) unique(x$PSM[x$isOutlier]))

message("\nSummary of PSM outliers:")
x <- sapply(psm_outliers,length)
data.table(Mixture=names(x), nOutliers=x) %>% knitr::kable()

# psm data for each mix in a list
psm_list <- msstats_psm %>% group_by(Mixture) %>% group_split()
names(psm_list) <- sapply(psm_list,function(x) unique(as.character(x$Mixture)))

# loop to remove outlier psm from each mixture
# NOTE: PSM with incomplete observations within a mixture (n=2) are removed
filt_list <- list()
for (mixture in names(psm_list)) {
	filt_psm <- psm_list[[mixture]] %>% 
		filter(PSM %notin% psm_outliers[[mixture]])
	filt_list[[mixture]] <- filt_psm
}


# collect results--outliers removed from each mixture
filt_psm <- do.call(rbind,filt_list)

stopifnot(!any(is.na(filt_psm$Intensity)))


## [2] summarize protein level data ----------------------------------------------
# Perform protein summarization for each run with MSstatsTMT.

# NOTE: my fork allows you to pass additional args to underlying MSstats 
# dataProcess function  -- speed things up by specifying the number of cores to 
# be used for parallel processing.

message("\nPerforming normalization and protein summarization using MSstatsTMT.")


n_cores <- parallel::detectCores() - 1

t0 <- Sys.time()

suppressMessages({ # verbosity
  msstats_prot <- proteinSummarization(filt_psm,
    method = "msstats", 
    remove_norm_channel = remove_norm_channel,
    global_norm = global_norm, # perform global norm using 'norm' condition
    MBimpute = MBimpute, 
    reference_norm = reference_norm,
    clusters = n_cores
  )
})


# This takes about 11 minutes for 8.5 k proteins with 23 cores
# FIXME: fix warnings messages about closing clusters.
proteins <- unique(as.character(msstats_prot$Protein))
message(
  "\nTime to summarize ", length(proteins), " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## [3] perform statistical comparisons with MSstatsTMT --------------------------

# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# fx <- formula(Abundance ~ 1 + Condition + (1|Mixture)) # lmerTest::lmer

# We specify Condition as Genotype.BioFraction for all intra-fraction
# comparisons. T-statistics are moderated using ebayes methods in limma.

# We perform the two analyses seperately so we can specify moderated = TRUE for 
# intra-BioFraction comparisons and moderated = FALSE for overall
# 'Mutant-Control' comparison.

message("\nAssessing protein-level comparisons with MSstatsTMT.")

t0 <- Sys.time()

## 1. 'intra-BioFraction' comparisons
suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results1 <- groupComparisonTMT(msstats_prot,
      msstats_contrasts,
      moderated = TRUE
    )
  })
})

# examine the results
results1 %>% group_by(Label) %>% 
	summarize(nSig = sum(adj.pvalue < FDR_alpha),.groups="drop") %>% 
	knitr::kable()

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 8 'intra-BioFraction' comparisons for ", length(proteins),
  " proteins: ", round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


## 2. 'Mutant-Control' comparison

t0 <- Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results2 <- groupComparisonTMT(msstats_prot,
      mut_vs_control,
      moderated = FALSE
    )
  })
})

# Adjust pvalues for multiple comparisons with Bonferroni method
results2$adj.pvalue <- p.adjust(results2$pvalue,method="bonferroni")

# examine the results
results2 %>% group_by(Label) %>% 
	summarize(nSig = sum(adj.pvalue < FDR_alpha),.groups="drop") %>% 
	knitr::kable()

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 'Mutant-Control' comparison for ", length(proteins),
  " proteins: ", round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")

msstats_results <- rbind(results1,results2)

## save results ---------------------------------------------------------------

# save msstats_prot -- the normalized protein data
myfile <- file.path(root, "data", "msstats_prot.rda")
save(msstats_prot, file = myfile, version = 2)

message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# save msstats_results
myfile <- file.path(root, "data", "msstats_results.rda")
save(msstats_results, file = myfile, version = 2)

message("\nSaved ", basename(myfile), " in ", dirname(myfile))
