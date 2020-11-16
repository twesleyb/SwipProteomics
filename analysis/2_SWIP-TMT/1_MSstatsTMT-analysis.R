#!/usr/bin/env Rscript 

# title: MSstatsTMT
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Options
FDR_alpha <- 0.05 
save_rda <- TRUE

## MSstatsTMT Options
moderated <- TRUE
MBimpute <- TRUE
rm_single <- TRUE
global_norm <- TRUE
reference_norm <- TRUE
remove_norm_channel <- TRUE


## prepare the working environment ---------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# load functions in root/R and make data in root/data accessible
# library(SwipProteomics)
devtools::load_all(root)

# load data in root/data
data(pd_psm) 
data(gene_map)
data(pd_annotation)
data(mutant_vs_control) # 'Mutant-Control' comparison
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


## Functions ------------------------------------------------------------------

detectOutliers <- function(mixture, psm_df, nbins=5, nSD=4, nComplete=2) {
  # this function identifies QC PSM-level outliers
  # NOTE: PSM with incomplete (n != nComplete) features are removed
  # We cannot assess reproducibility of a PSM if it was only quantified once
  # collect psm with incomplete (missing) data
  is_incomplete <- psm_df %>% dplyr::filter(Mixture == mixture & Condition == "Norm") %>%
	  group_by(PSM) %>% mutate(nObs=length(PSM)) %>% ungroup() %>% filter(nObs != nComplete) %>%
	  select(PSM) %>% unlist() %>% unique()
  # subset, keep psm with complete features
  filt_df <- psm_df %>% dplyr::filter(Mixture == mixture & Condition == "Norm") %>%
	  group_by(PSM) %>% mutate(nObs=length(PSM)) %>% ungroup()  %>%
          # drop PSM with missing features
	  filter(nObs == nComplete) %>% 
	  # calculate mean of log2 Intensity for binning psm
	  group_by(Mixture,PSM) %>% mutate(meanQC=mean(log2(Intensity))) %>% 
	  ungroup()
  # group ratio data into intensity bins
  breaks <- stats::quantile(filt_df$meanQC, seq(0, 1, length.out=nbins+1),
  		   names=FALSE,na.rm=TRUE)
  filt_df <- filt_df %>% 
  	mutate(bin = cut(meanQC, breaks, labels=FALSE, include.lowest=TRUE))
  # compute log ratio (difference) of QC measurments
  filt_df <- filt_df %>% group_by(Mixture,PSM) %>% 
  	mutate(ratioQC = diff(log2(Intensity))) %>% ungroup()
  # summarize mean and SD of each bin
  filt_df <- filt_df %>% group_by(bin) %>% 
  	mutate(meanRatio = mean(ratioQC,na.rm=TRUE),
  	       binSD = sd(ratioQC), .groups="drop") %>% ungroup()
  # flag outliers
  filt_df <- filt_df %>% mutate(isLow = ratioQC < meanRatio - nSD*binSD)
  filt_df <- filt_df %>% mutate(isHigh = ratioQC > meanRatio + nSD*binSD)
  filt_df <- filt_df %>% mutate(isOutlier = isLow | isHigh)
  # summary
  summary_df <- filt_df %>% group_by(bin) %>%
  	summarize(meanIntensity = mean(Intensity),
  		  meanRatio = mean(ratioQC,na.rm=TRUE),
  		  binSD = sd(ratioQC), 
  		  nOut = sum(isOutlier), .groups="drop")
  # status
  if (length(is_incomplete) > 0) {
    warning(length(is_incomplete),
	    " PSM with incomplete features in mixture ", mixture,
	    " were removed.", call.=FALSE)
  }
  out_df <- filt_df %>% select(Mixture,PSM, bin, meanQC, ratioQC, binSD, isOutlier)
  return(out_df %>% filter(isOutlier))
} # EOF


## [1] convert to msstats format -----------------------------------------------
# Proteins with a single feature (i.e. peptide) are removed.

message("\nConverting PD PSM-level data into MSstatsTMT format.")

# aprox 7 minutes for all proteins
suppressMessages({ # verbosity
  msstats_psm <- PDtoMSstatsTMTFormat(pd_psm,
    pd_annotation, # see 0_PD-data-preprocess.R
    which.proteinid = "Master.Protein.Accessions",
    rmProtein_with1Feature = rm_single
  )
})


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
outlier_list <- lapply(mix, detectOutliers, msstats_psm, nbins=5, nSD=4)
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
# Perform protein summarization for each run.

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
    global_norm = global_norm, # perform global normalization using 'norm'
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

# 1. 'intra-BioFraction' comparisons
suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results1 <- groupComparisonTMT(msstats_prot,
      msstats_contrasts,
      moderated = moderated
    )
  })
})


results1 %>% subset(Protein %in% washc_prots)

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 8 'intra-BioFraction' comparisons for ", length(proteins),
  " proteins: ", round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


t0 <- Sys.time()

# 1. 'Mutant-Control' comparison
suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results2 <- groupComparisonTMT(msstats_prot,
      mut_vs_control,
      moderated = FALSE
    )
  })
})

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 'Mutant-Control' comparison for ", length(proteins),
  " proteins: ", round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")



## clean-up and combine results -----------------------------------------------

# combine results
msstats_results <- rbind(results1,results2)


## save results ---------------------------------------------------------------

if (save_rda) {

  ## save impute_prot -- the normalized, imputed protein data
  #myfile <- file.path(root, "data", "impute_prot.rda")
  #save(impute_prot, file = myfile, version = 2)
  #message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # save msstats_prot -- the normalized protein data
  myfile <- file.path(root, "data", "msstats_prot.rda")
  save(msstats_prot, file = myfile, version = 2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # save results
  myfile <- file.path(root, "data", "msstats_results.rda")
  save(msstats_results, file = myfile, version = 2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # save sigprots
  myfile <- file.path(root,"data","sigprots.rda")
  save(sigprots,file=myfile,version=2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

}
