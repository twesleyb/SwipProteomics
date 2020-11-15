#!/usr/bin/env Rscript 

# title: MSstatsTMT
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root <- "~/projects/SwipProteomics"

## Options
FDR_alpha <- 0.05 
moderated <- TRUE
save_rda <- TRUE
MBimpute <- TRUE
rm_single <- TRUE
global_norm <- TRUE
reference_norm <- TRUE

# rm proteins with a single feature


## prepare the working environment ---------------------------------------------

# load renv
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  # my MSstats forks:
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})


## NOTE: my fork attempts to remove much of MSstats's verbosity.
## NOTE: you may find some debug.log files on your computer still, their source
## has not been found yet.
## twesleyb/MSstatsTMT includes access to the internal functions used by
## MSstatsTMT to fit protein-wise models and perform statistical comparisons.
## MSstatsTMT is a wrapper around MSstats. My fork allows you to pass arguments
## for parallel processing to MSstats::proteinSummarization to speed things up.


## load data ------------------------------------------------------------------

# load functions in root/R and make data in root/data accessible
# library(SwipProteomics)
devtools::load_all()

# load data in root/data
data(pd_psm) # 46 mb
data(gene_map)
data(pd_annotation)
data(msstats_contrasts)

# NOTE: msstats_contrasts is a matrix specifying pairwise contrasts between all
# 'BioFraction.Control' and 'BioFraction.Mutant' Conditions.


## Functions ------------------------------------------------------------------

cleanProt <- function(df) {
  # clean-up msstats_prot results data.frame
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
  # clean-up msstat_results data.frame
  idx <- match(df$Protein, gene_map$uniprot)
  df$Symbol <- gene_map$symbol[idx]
  df$Entrez <- gene_map$entrez[idx]
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


## 0. Create a contrast -------------------------------------------------------

# add a contrast vector specifying the mutant-control comparison to
# msstats_contrasts matrix (see 0_*.R)

mut_vs_control <- matrix(c(
  -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7,
  1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7
), nrow = 1)
row.names(mut_vs_control) <- "Mutant-Control"
colnames(mut_vs_control) <- colnames(msstats_contrasts)

msstats_contrasts <- rbind(msstats_contrasts, mut_vs_control)


## [1] convert to msstats format -----------------------------------------------
# Proteins with a single feature (i.e. peptide) are removed.

message("\nConverting PD PSM-level data into MSstatsTMT format.")

t0 <- Sys.time()

suppressMessages({ # verbosity
  msstats_psm <- PDtoMSstatsTMTFormat(pd_psm,
    pd_annotation, # see 0_PD-data-preprocess.R
    which.proteinid = "Master.Protein.Accessions",
    rmProtein_with1Feature = rm_single
  )
})

# aprox 7 minutes for all proteins
message(
  "\nTime to pre-process ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## [added] QC-based PSM filtering ------------------------------------------------

# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated.
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, measuremnts that are outside
# +/- nSD x standard deviations from the bin's mean are removed.

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
filt_msstats_psm <- do.call(rbind,filt_list)


## [2] summarize protein level data ----------------------------------------------
# Perform protein summarization for each run.

# NOTE: my fork allows you to pass additional args to underlying MSstats 
# dataProcess function  -- speed things up by specifying the number of cores to 
# be used for parallel processing.

n_cores <- parallel::detectCores() - 1

message("\nPerforming normalization and protein summarization using MSstatsTMT.")

t0 <- Sys.time()

suppressMessages({ # verbosity
  msstats_prot <- proteinSummarization(filt_msstats_psm,
    method = "msstats",
    remove_norm_channel = TRUE,
    global_norm = global_norm, # perform global normalization using 'norm'
    MBimpute = MBimpute, 
    reference_norm = reference_norm,
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
# fx <- formula(Abundance ~ 1 + Condition + (1|Mixture)) # lmerTest::lmer

# We specify Condition as Genotype.BioFraction for all intra-fraction
# comparisons. T-statistics are moderated using ebayes methods in limma.

message("\nAssessing protein-level comparisons with MSstatsTMT.")

t0 <- Sys.time()

suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    msstats_results <- groupComparisonTMT(msstats_prot,
      msstats_contrasts,
      moderated = TRUE
    )
  })
})

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 9 pairwise comparisons for ", nprot, " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


### [added] impute protein-level missing values ---------------------------------
#
## There are multiple levels of missingness:
## * missingness within a Run -- handled by MSstats during protein processing
## * missingness between mixtures
#
## there are proteins with only a couple missing values that we have to discard
## because missing values are not tolerated in a number of the subsequent steps
## MSstats should have imputed these (I thought), but here we do so with KNN for
## MNAR data.
#
## FIXME: need to visualize this
#message("\nImputing protein-level missing values with KNN.")
#
## cast the data into a matrix
#dm <- msstats_prot %>% 
#  reshape2::dcast(Protein ~ Mixture + Channel + Condition,
#		value.var="Abundance") %>%
#  as.data.table() %>% as.matrix(rownames="Protein")
#
#
## we will impute up to 50% missingness--a protein must be identified in at least
## 1 mixture/experiment.
#idx <- apply(dm,1,function(x) sum(is.na(x)) > 0.5 * ncol(dm))
#if (sum(idx)>0) {
#  warning(sum(idx)," rows (proteins) with > 50% missingness cannot be imputed.")
#}
#
## impute, suppress output from cat with sink()
#sink(tempfile())
#dm_knn <- impute::impute.knn(dm[!idx,],k=10,colmax=0.8,rowmax=0.5)$data
#sink()
#
## there should be no missing values
#stopifnot(!any(is.na(dm_knn)))
#
## tidy the imputed data
#df_knn <- reshape2::melt(dm_knn)
#colnames(df_knn) <- c("Protein","Sample","Abundance")
#df_knn$Mixture <- sapply(strsplit(as.character(df_knn$Sample),"_"),"[",1)
#df_knn$Channel <- sapply(strsplit(as.character(df_knn$Sample),"_"),"[",2)
#df_knn$Condition <-sapply(strsplit(as.character(df_knn$Sample),"_"),"[",3)
#
#df_knn$censor <- TRUE # censor these values when computing correlation stats
#
## merge with msstats_prot, write over Abundance, keep proteins with complete obs
#cols <- intersect(colnames(msstats_prot),colnames(df_knn))
## NOTE: use right_join here!
#impute_prot <- msstats_prot %>% 
#	filter(Protein %in% df_knn$Protein) %>% 
#	right_join(df_knn,by=cols)
#
#
## proceed with imputed data -- there should be no missing observations
## proteins with missingness have been imputed or removed
## status
#nprots <- length(unique(impute_prot$Protein))
#message("Final number of proteins with complete observations: ", 
#	formatC(nprots,big.mark=","))
#
#stopifnot(!(any(is.na(impute_prot$Abundance))))


## clean-up and combine results -----------------------------------------------

# clean-up protein results
msstats_prot <- cleanProt(msstats_prot)

# collect a list of results 
results_list <- cleanResults(msstats_results) %>% as.data.table() %>%
	group_by(Contrast) %>% group_split()
class(results_list) <- "list"

# sort list
namen <- sapply(results_list,function(x) unique(x$Contrast))
names(results_list) <- gsub("Mutant.F[0-9]{1,2}-Control\\.","",namen)
idx <- c("F4", "F5", "F6", "F7", "F8", "F9", "F10","Mutant-Control")
results_list <- results_list[idx]


# combine into list to be saved as excel
results_list <- c(
  "Normalized Protein" = list(msstats_prot),
  results_list
)

# save results as excel document
myfile <- file.path(root, "tables", "S2_SWIP_TMT_Results.xlsx")
write_excel(results_list, myfile)


## summarize significant results ----------------------------------------------

message("\nTotal instances of significant change: ",
	sum(sapply(results_list,function(x) sum(x$FDR<FDR_alpha))))

message("\nSummary of significant proteins for each contrast:")
sapply(results_list[-1],function(x) sum(x$FDR<FDR_alpha)) %>% 
	t() %>% knitr::kable()

# proteins with sig change
sigprots <- msstats_results %>% filter(adj.pvalue<FDR_alpha) %>% 
	select(Protein) %>% unlist() %>% as.character() %>% unique()


## save results ---------------------------------------------------------------

if (save_rda) {

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
