#!/usr/bin/env Rscript 
# title: MSstatsTMT
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root <- "~/projects/SwipProteomics"

## Options
nprot <- "all" # the number of proteins to be analyzed or 'all'
FDR_alpha <- 0.05 # FDR threshold for significance
save_rda <- TRUE


## prepare the working environment ---------------------------------------------

# load renv
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})


## NOTE: my fork attempts to remove much of MSstats's verbosity as well as
## includes access to the internal functions used by MSstatsTMT to fit
## protein-wise models and perform statistical comparisons between groups.
## MSstatsTMT is a wrapper around MSstats. My fork allows you to pass
## arguments for parallel processing to proteinSummarization to speed things up.


## load data ------------------------------------------------------------------

# load functions in root/R and make data in root/data accessible
suppressWarnings({
  devtools::load_all()
})

# load data in root/data
data(pd_psm) # 46 mb
data(gene_map)
data(pd_annotation)
data(msstats_contrasts)
# NOTE: msstats_contrasts is a matrix indicating pairwise contrasts between all
# BioFraction.Control and BioFraction.Mutant


## Functions ------------------------------------------------------------------

cleanProt <- function(df) {
  # clean-up msstats_prot
  # add BioFraction, Genotype, and Subject columns
  idx <- match(df$Protein, gene_map$uniprot)
  Symbol <- gene_map$symbol[idx]
  Entrez <- gene_map$entrez[idx]
  df <- tibble::add_column(df, Symbol, .after = "Protein")
  df <- tibble::add_column(df, Entrez, .after = "Symbol")
  geno <- sapply(strsplit(as.character(df$Condition), "\\."), "[", 1)
  biof <- sapply(strsplit(as.character(df$Condition), "\\."), "[", 2)
  subject <- as.numeric(interaction(df$Mixture, geno))
  df$Genotype <- factor(geno, levels = c("Mutant", "Control"))
  df$BioFraction <- factor(biof,
    levels = c("F4", "F5", "F6", "F7", "F8", "F9", "F10")
  )
  df$Subject <- as.factor(subject)
  df <- df %>% select(
    "Run", "Mixture", "TechRepMixture", "Channel", "Condition",
    "BioFraction", "Genotype", "Subject", "BioReplicate",
    "Protein", "Symbol", "Entrez", "Abundance"
  )
  return(df)
}


cleanResults <- function(df) {
  # clean-up msstat_results
  idx <- match(df$Protein, gene_map$uniprot)
  df$Symbol <- gene_map$symbol[idx]
  df$Entrez <- gene_map$entrez[idx]
  df <- df %>% filter(is.na(issue))
  df$issue <- NULL
  df <- df %>%
    select(
      Label, Protein, Entrez,
      Symbol, log2FC, pvalue, adj.pvalue, SE, DF
    ) %>%
    unique() %>%
    arrange(pvalue)
  colnames(df)[colnames(df)=="Label"] <- "Contrast"
  colnames(df)[colnames(df)=="pvalue"] <- "Pvalue"
  colnames(df)[colnames(df)=="adj.pvalue"] <- "FDR"
  return(df)
}


detectOutliers <- function(mixture, psm_df, nbins=5, nSD=4, nComplete=2) {
  psm_df <- psm_df %>% dplyr::filter(Mixture == mixture & Condition == "Norm") %>%
	  group_by(PSM) %>% mutate(nObs=length(PSM)) %>% ungroup() %>%
          # drop PSM with missing features
	  filter(nObs == nComplete) %>% 
	  # calculate mean of log2 Intensity for binning psm
	  group_by(Mixture,PSM) %>% mutate(meanQC=mean(log2(Intensity))) %>% 
	  ungroup()
  # Group ratio data into intensity bins.
  breaks <- stats::quantile(psm_df$meanQC, seq(0, 1, length.out=nbins+1),
  		   names=FALSE,na.rm=TRUE)
  psm_df <- psm_df %>% 
  	mutate(bin = cut(meanQC, breaks, labels=FALSE, include.lowest=TRUE))
  # compute log ratio (difference) of QC measurments
  psm_df <- psm_df %>% group_by(Mixture,PSM) %>% 
  	mutate(ratioQC = diff(log2(Intensity))) %>% ungroup()
  # summarize mean and SD of each bin
  psm_df <- psm_df %>% group_by(bin) %>% 
  	mutate(meanRatio = mean(ratioQC,na.rm=TRUE),
  	       binSD = sd(ratioQC), .groups="drop") %>% ungroup()
  # flag outliers
  psm_df <- psm_df %>% mutate(isLow = ratioQC < meanRatio - nSD*binSD)
  psm_df <- psm_df %>% mutate(isHigh = ratioQC > meanRatio + nSD*binSD)
  psm_df <- psm_df %>% mutate(isOutlier = isLow | isHigh)
  # summary
  summary_df <- psm_df %>% group_by(bin) %>%
  	summarize(meanIntensity = mean(Intensity),
  		  meanRatio = mean(ratioQC,na.rm=TRUE),
  		  binSD = sd(ratioQC), 
  		  nOut = sum(isOutlier), .groups="drop")
  out_df <- psm_df %>% select(Mixture,PSM, bin, meanQC, ratioQC, binSD, isOutlier)
  return(out_df %>% filter(isOutlier))
} # EOF


## [0] subset the data ------------------------------------------------------------
# NOTE: we have previously made sure Master.Protein.Accessions is a column of
# (single) Uniprot Identifiers. NOTE: Only run once!

# all proteins
proteins <- unique(as.character(pd_psm$Master.Protein.Accessions))

# define proteins to be analyzed:
if (nprot > length(proteins) | nprot == "all") {
  # using all proteins
  message(
    "\nAnalyzing 'all' (N=",
    formatC(length(proteins), big.mark = ","), ") proteins."
  )
  msstats_input <- pd_psm
} else {
  # subset the data
  message("\nAnalyzing a subset of the data (n=", nprot, ") proteins.")
  prots <- sample(proteins, nprot)
  pd_filt <- pd_psm %>%
    filter(Master.Protein.Accessions %in% prots)
  msstats_input <- pd_filt
} # EIS


## [1] convert to msstats format -----------------------------------------------
# Proteins with a single feature (i.e. peptide) are removed.

message("\nConverting PD PSM-level data into MSstatsTMT format.")

t0 <- Sys.time()

suppressMessages({ # verbosity
  msstats_psm <- PDtoMSstatsTMTFormat(msstats_input,
    pd_annotation,
    which.proteinid = "Master.Protein.Accessions",
    rmProtein_with1Feature = TRUE
  )
})

# aprox 7 minutes for all proteins
message(
  "\nTime to pre-process ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)

## [added] QC-based PSM filtering ----------------------------------------------------

# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated.
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, measuremnts that are outside
# +/- nSD x standard deviations from the bin's mean are removed.

# This is essential because our normalization strategy depends upon the SPQC
# samples. We cannot perform normalization to SPQC PSM that are highly variable.

# collect all psm outliers for each Mixture
mix <- c("M1","M2","M3")
outlier_list <- lapply(mix, detectOutliers, msstats_psm, nbins=5, nSD=4)
names(outlier_list) <- mix

# all psm outliers
psm_outliers <- sapply(outlier_list,function(x) unique(x$PSM[x$isOutlier]))

sapply(psm_outliers,length)

# loop to remove outlier psm from each mixture
filt_list <- list()

psm_list <- msstats_psm %>% group_by(Mixture) %>% group_split()
names(psm_list) <- sapply(psm_list,function(x) unique(as.character(x$Mixture)))

for (mixture in names(psm_list)) {
	filt_psm <- psm_list[[mixture]] %>% filter(PSM %notin% psm_outliers[[mixture]])
	filt_list[[mixture]] <- filt_psm
}


# collect results--outliers removed from each mixture
filt_msstats_psm <- do.call(rbind,filt_list)


## [2] summarize protein level data ----------------------------------------------
# Perform protein summarization for each run.

# NOTE: my fork allows you to pass additional args to underlying MSstats 
# dataProcess function  -- speed things up by specifying the number of cores to 
# be used for parallel processing.

n_cores <- parallel::detectCores() - 1

message("\nPerforming normalization and protein summarization using MSstatsTMT.")

t0 <- Sys.time()

# NOTE: Don't remove normalization channel!
# MSstats does not perform any batch (aka Mixture or Experiment) normalization.
# The variance due to random variability of Mixture is dealt with by the
# mixed-models used to assess protein differential abundance. However,
# we wish to contruct a covariation matrix from proteins that covary together in
# subcellular space. We should remove any variability that can be explained by
# Mixture effects before we do this.

suppressMessages({ # verbosity
  msstats_prot <- proteinSummarization(filt_msstats_psm, #  or msstats_psm
    method = "msstats",
    remove_norm_channel = TRUE,
    global_norm = TRUE,
    MBimpute = TRUE,
    reference_norm = TRUE,
    clusters = n_cores
  )
})

# This takes about 11 minutes for 8.5 k proteins with 23 cores
# FIXME: fix warnings messages about closing clusters.
message(
  "\nTime to summarize ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## [3] perform intrafraction statistical comparisons --------------------------

# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# lmerTest::lmer(Abundance ~ (1|Mixture) + Condition)
fx0 <- formula(Abundance ~ 0 + Condition + (1|Mixture))

# We specify Condition as Genotype.BioFraction for all intra-fraction
# comparisons. T-statistics are moderated using ebayes methods in limma.

message("\nAssessing intrafraction comparisons with MSstatsTMT.")

t0 <- Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity
    msstats_results <- groupComparisonTMT(msstats_prot,
      msstats_contrasts,
      moderated = TRUE
    )
  })
})

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform intra-Biofraction comparisons for ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## [4] perform statistical comparisons for 'Control-Mutant' comparison ---------

# create a contrast for assessing difference between 'Control' and 'Mutant' coef:
alt_contrast <- matrix(c(
  -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7,
  1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7
), nrow = 1)
row.names(alt_contrast) <- "Mutant-Control"
colnames(alt_contrast) <- levels(msstats_prot$Condition)

message("\nAssessing 'Mutant-Control' comparison with MSstatsTMT.")

t0 <- Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity
    res2 <- MSstatsTMT::groupComparisonTMT(
      data = msstats_prot,
      contrast.matrix = alt_contrast,
      moderated = FALSE # no moderation
    )
  })
})

# This takes about ?? minutes for 8.5 k proteins
message(
  "\nTime to perform comparison for ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## clean-up and combine results -----------------------------------------------

# clean-up protein results
msstats_prot <- cleanProt(msstats_prot)

# combine results
tmp_list <- cleanResults(msstats_results) %>%
  group_by(Contrast) %>%
  group_split()
names(tmp_list) <- sapply(tmp_list, function(x) unique(x$Contrast))

# simplify names
names(tmp_list) <- sapply(strsplit(names(tmp_list), "\\."), "[", 3)

# sort list
tmp_list <- tmp_list[c("F4", "F5", "F6", "F7", "F8", "F9", "F10")]

# fix column names
tmp_list <- lapply(tmp_list,function(x) {
  colnames(x)[colnames(x) == "pvalue"] <- "Pvalue"
  colnames(x)[colnames(x) == "adj.pvalue"] <- "FDR"
  return(x)
})


# sort msstats_prot columns 
tmp_df <- msstats_prot %>% select(Mixture,Channel,Genotype,BioFraction,
				  Condition,Subject,Protein,Symbol,Entrez,Abundance)

# combine into list to be saved as excel
results_list <- c(
  "Normalized Protein" = list(tmp_df),
  tmp_list, # intrafraction results
  "Mutant-Control" = list(cleanResults(res2))
)


# save results as excel document
myfile <- file.path(root, "tables", "S2_SWIP_TMT_Results.xlsx")
write_excel(results_list, myfile)


## summarize significant results ----------------------------------------------

message("\nTotal instances of significant change: ",
	sum(sapply(tmp_list,function(x) sum(x$FDR<FDR_alpha))))

message("\nSummary of significant proteins:")
sapply(tmp_list,function(x) sum(x$FDR<FDR_alpha)) %>% t() %>% knitr::kable()

results_list[["Mutant-Control"]] %>% 
	summarize(Contrast = unique(Contrast), nSig = sum(FDR<FDR_alpha)) %>% 
	knitr::kable()


## save results ---------------------------------------------------------------

if (save_rda) {

  # save the model formula
  myfile <- file.path(root, "data", "fx0.rda")
  save(fx0,file=myfile,version=2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # save msstats_prot -- the normalized protein data
  myfile <- file.path(root, "data", "msstats_prot.rda")
  save(msstats_prot, file = myfile, version = 2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # save results
  msstats_results <- cleanResults(rbind(msstats_results, res2))
  myfile <- file.path(root, "data", "msstats_results.rda")
  save(msstats_results, file = myfile, version = 2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

  # quick, save sigprots
  sigprots <- msstats_results %>% 
	filter(Contrast == 'Mutant-Control' & FDR < 0.05) %>%
	select(Protein) %>% unlist() %>% as.character() %>% unique()
  myfile <- file.path(root,"data","sigprots.rda")
  save(sigprots,file=myfile,version=2)
  message("\nSaved ", basename(myfile), " in ", dirname(myfile))

}
