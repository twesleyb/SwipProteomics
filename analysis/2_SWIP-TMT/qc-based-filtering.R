#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: Preprocessing and statistical analysis of Swip TMT proteomics.
#' authors: Tyler W A Bradshaw
#' ---

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.

## INPUT data is in root/data/TMT_.zip.
# The zipped directory contains sample meta data and raw peptide data from PD:

# NOTE: this data is no longer being used as input for the rest of the analysis
# but we consider it as an example. 
zip_file = "TMT.zip"
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

## OPTIONS:
fig_height = 5.0 # Default height of figures (in).
fig_width = 5.0 # Default width of figures (in).

oldham_threshold = 2.5 # Sample connectivity threshold for detecting outliers.
# NOTE: No sample outliers were detected.

FDR_alpha = 0.1 # FDR threshold for differential abundance.
BF_alpha = 0.05 # Bonferroni threshold for differential abundance.
logFC_threshold = c(lwr=log2(0.8), upr=log2(1.2)) # FC threshold.

set.seed = as.numeric(Sys.time()) # seed for random operations.

## OUTPUT saved in root/.../tables:
# * Swip_TMT_Results.xlsx - an excel spreadsheet that contains:
#     - Sample meta data.
#     - Gene/protein identifiers.
#     - Raw peptide data.
#     - Final normalized data (log2 transformed).
#     - Statistical results for intrafraction contrasts (F4-10).
#     - Statistical results for WT v MUT contrast.

## OUTPUT for R package in root/data.
# Key datasets are saved in the data directory as read-only objects.
# * tmt_protein.rda
# NOTE: can be loaded in R with data(tmt_protein).

## OUTPUT for downstream analysis in root/rdata/
# Input files/temporary files that are passed to other scripts including
# the Python Leidenalg clustering script are saved in root/rdata.

## Order of data processing operations:
# * Load the raw peptide intensity data from PD 2.2.
# * Initial sample loading normalization.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility - remove QC outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization.
# * IRS normalization - equalize protein measurements from different peptides.
# * Sample pool normalization - assumes that mean of pooled QC technical
#   replicates is equal to mean of all biological replicates.
# * Protein level filtering - remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not
#   reproducible (exhibit high inter-experimental variablility across the 3x
#   biological replicates).
# * Assess intra-fraction differential abundance with a general linear model.
#   GLM: [Protein Abundance] ~ 0 + Fraction.Genotype.
# * Assess WT v MUT differential abundance across all fractions with a
#   GLM: [Protein Abundance] ~ Fraction + Genotype.

## NOTE:
# The dependencies for this project were managed within a conda virtual
# environment into which R and renv, an R package for managing R libraries,
# were installed. All additional R dependencies were installed using
# renv. The anRichment package was installed from source (see installation
# script in root/).

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) {
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE)))
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) {
		here <- dirname(here)
	}
	root <- here
	return(root)
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis.

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE)

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(data.table) # For working with tables.
	library(getPPIs) # For PPIs and mapping gene identifiers.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(rootdir, "data") # Key pieces of data saved as rda
rdatdir <- file.path(rootdir, "rdata") # Temporary data files
downdir <- file.path(rootdir, "downloads") # Misc downloads/temporary files

# Create output directories if necessary.
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

#---------------------------------------------------------------------
## Load the raw data and sample info.
#---------------------------------------------------------------------

# Extract the raw TMT data from zipped file.
myfile <- file.path(datadir,zip_file)
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# Read TMT.csv data into R with data.table::fread.
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_data)
peptides <- fread(myfile)

# Load sample information.
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- fread(myfile)

# Format Cfg force column -- this is the Force in g's used to obtain
# the subcellular fraction.
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")

rm(list="myfile") # remove temporary variables.

#---------------------------------------------------------------------
## Map all Uniprot accession numbers to stable entrez IDs.
#---------------------------------------------------------------------

# Map Uniprot IDs to Entrez using MGI batch query function.
# Map Entrez IDs to gene symbols using getPPIs::getIDs().
message("\nCreating gene identifier map.")

# First, remove any non-mouse proteins from the data.
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# Remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# Collect all Uniprot IDs.
uniprot <- unique(peptides$Accession)

# Map Uniprot IDs to Entrez using online MGI batch query function.
# This takes a couple minutes because the function currently
# downloads the MGI data each time the function is called.
entrez <- mgi_batch_query(uniprot,quiet=TRUE)
names(entrez) <- uniprot

# Map any remaining missing IDs by hand.
message("Mapping missing IDs by hand.\n")
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(P05214=22144, P0CG14=214987, P84244=15078)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }

# Map entrez ids to gene symbols using twesleyb/getPPIs.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Create gene identifier mapping data.table.
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")

rm(list=c("ig_prots","uniprot","entrez","missing",
	  "mapped_by_hand", "check","symbols"))

#---------------------------------------------------------------------
## Tidy-up the input data from Proteome Discover.
#---------------------------------------------------------------------

# Convert PD df into tidy df.
# Samples should contain the following columns:
# Treatment, Channel, Sample, Experiment
message("\nLoading raw data from Proteome Discover (PD2.2).")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,id.vars=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

# Summary of peptide/protein quantification:
message(paste("\nSummary of initial peptide/protein quantification",
	      "after removing contaiminants:"))
n_samples <- length(unique(tidy_peptide$Sample))
n_proteins <- length(unique(tidy_peptide$Accession))
n_peptides <- length(unique(tidy_peptide$Sequence))

df <- data.frame("Samples"=as.character(n_samples),
		 "Proteins"=formatC(n_proteins,big.mark=","),
		 "Peptides"=formatC(n_peptides,big.mark=","))
knitr::kable(df)

rm(list=c("cols","df","n_samples","n_proteins","n_peptides"))

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

# Perform sample normalization. Normalization is done for each
# experiment independently (group by Experiment:Sample).
# NOTE: Grouping by Experiment, Channel does't work because
# Sample TMT channels (e.g. 126N were used in different expirements).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))

# Check the data before SL:
#message("\nPeptide data before sample loading normalization:")
#check_SL(tidy_peptide)

# Check the data after SL:
#essage("\nPeptide data after sample loading normalization:")
#heck_SL(sl_peptide) # Equal within an experiment.

#---------------------------------------------------------------------
## Impute missing peptide values.
#---------------------------------------------------------------------

# Impute missing peptide values with k-nearest neighbors (KNN) algorithm.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
# NOTE: KNN imputing is done for each experimental group seperately.
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC",quiet=FALSE)

#---------------------------------------------------------------------
## Examine reproducibility of QC measurements.
#---------------------------------------------------------------------
# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated.
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, measuremnts that are outside
# +/- 4x standard deviations from the bin's mean are removed.

message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide,controls="SPQC",quiet=FALSE)

# unpack 
tp <- imputed_peptide
controls = "SPQC"
nbins = 5
nSD = 4


filtQC <- function(tp,controls="SPQC",nbins=5,nSD=4,quiet=TRUE){
	# Remove peptides with highly variable QC measurments.
	# Calculate ratios of QC peptides, grouped by Experiment.

	# Imports.
	suppressPackageStartupMessages({
		library(dplyr)
		library(data.table)
	})

	ratio_data <- tp %>% 
		# get the control data
		filter(Treatment == controls) %>% 
		# group by PSM within an experiment
		group_by(Experiment,Accession,Sequence,Modifications) %>% 
		# calculate ratio of QC samples
		dplyr::summarize(Ratio = log2(Intensity[2]) - log2(Intensity[1]),
			  Mean = log2(mean(Intensity)),
			  N = length(Intensity),
			  nMissing = sum(is.na(Intensity)),
			  Remove = any(is.na(Intensity)), .groups="drop")

	# remove QC with any missing intensity 
	ratio_data <- ratio_data %>% filter(!Remove)

	# Group ratio data into intensity bins.
	breaks <- quantile(ratio_data$Mean,seq(0,1,length.out=nbins+1),
			   names=FALSE,na.rm=TRUE)
	ratio_data$Bin <- cut(ratio_data$Mean,breaks,
			      labels=FALSE,include.lowest=TRUE)

	# Summarize intensity bins.
	# this mostly for inspection
	ratio_df <- ratio_data %>% group_by(Bin) %>% 
		summarize("Median"= median(Mean),
			  "Mean"= mean(Ratio,na.rm=TRUE),
			  "Std" = sd(Ratio),
			  "N" = sum(!is.na(Ratio)),
			  "Min" = mean(Ratio,na.rm=TRUE)-(nSD*sd(Ratio)),
			  "Max" = mean(Ratio,na.rm=TRUE)+(nSD*sd(Ratio)),
			  .groups="drop") 

	# Determine if QC measurement is outside percision limits.
	ratio_data$Min <- ratio_df$Min[ratio_data$Bin]
	ratio_data$Max <- ratio_df$Max[ratio_data$Bin]
	out_low <- ratio_data$Ratio < ratio_data$Min 
	out_high <- ratio_data$Ratio > ratio_data$Max
	out <- out_low | out_high
	ratio_data$isOutlier <- out

	# Summarize number of outlies per bin.
	nOutliers <- ratio_data %>% group_by(Bin) %>% 
		summarize(n=sum(isOutlier),.groups="drop")
	ratio_df$nOutliers <- nOutliers[["n"]]

	# Collect outlier peptides.
	data_outliers <- ratio_data %>% filter(isOutlier)
	outlier_peptides <- paste(data_outliers$Experiment,
				  data_outliers$Accession,
				  data_outliers$Sequence,
				  data_outliers$Modifications,sep="_")

	# Remove outlier peptides from data.
	ids <- paste(tp$Experiment,tp$Accession,
		     tp$Sequence,tp$Modifications,sep="_")
	is_outlier <- ids %in% outlier_peptides
	tp_filt <- tp %>% filter(!is_outlier)
	tp_filt <- as.data.frame(tp_filt)

	# Status report.
	if (!quiet){
		total_out <- sum(out,na.rm=TRUE)
		message(paste("Total number of",controls,
			      "outlier peptides identified:",total_out))
	}
	# Return tidy data.
		return(tp_filt)
}


## do it with msstats_psm level data
root <- getrd()
myfile <- file.path(root,"rdata","msstats_psm.rda")
load(myfile)

colnames(msstats_psm)

x = msstats_psm %>% group_by(PSM) %>% group_split()
length(x) # 109134 
# total psms?

n = sample(length(x),1)
y = x[[n]]
#dim(y)[1] == 32 # not identified in run 3
# dim
# may be incomplete -- only identified in 1,2, or 3 of the mixtures...
# really we want the post-imputed data from msstats... but this only exists
# inside msstats
# we can consider only complete cases...

y  = sapply(x,nrow)
idx <- y == 48 # cuts things down to 1/3
mylist = x[which(idx)]
y$Condition # remember there are 32 RUNS for our 48 fractions???

df = msstats_psm %>% group_by(PSM) %>% mutate(nObs=length(PSM))
subdf = df %>% filter(nObs == 48)
mylist <- subdf %>% group_by(PSM) %>% group_split()

x = mylist[[1]]

# now proceed with outlier detection
# do all the work to identify outliers
tmp_df <- msstats_psm %>% filter(Condition == "Norm") %>% 
	group_by(Mixture, PSM) %>% mutate(nObs=length(PSM)) %>% ungroup()
psm_df <- tmp_df %>% filter(nObs == nComplete) %>% 
	group_by(Mixture,PSM) %>% mutate(meanQC=mean(log2(Intensity))) %>% 
	ungroup()

# split data into groups by mixture
psm_list <- psm_df %>% group_by(Mixture) %>% group_split()

# for each mixture, id outlier psm
psm_outliers <- lapply(psm_list, detectOutliers)
names(psm_outliers) <- sapply(psm_list,function(x) {
				      unique(as.character(x$Mixture))})

# save
myfile = file.path(root,"data","psm_outliers.rda")
save(psm_outliers,file=myfile,version=2)


detectOutliers <- function(psm_df) {
  nComplete <- 2 # techRep Norm
  nSD <- 4
  # this should be done for each experiment/mixture
  # Group ratio data into intensity bins.
  breaks <- stats::quantile(psm_df$meanQC, seq(0, 1, length.out=nbins+1),
  		   names=FALSE,na.rm=TRUE)
  psm_df <- psm_df %>% 
  	mutate(bin = cut(meanQC, breaks, labels=FALSE, include.lowest=TRUE))
  # compute ratio
  psm_df <- psm_df %>% group_by(Mixture,PSM) %>% 
  	mutate(ratioQC=log2(Intensity[1]) - log2(Intensity)[2]) %>% ungroup()
  # summarize bins
  psm_df <- psm_df %>% group_by(bin) %>% 
  	mutate(meanRatio = mean(ratioQC,na.rm=TRUE),
  	       sdRatio = sd(ratioQC), .groups="drop")  %>% ungroup()
  # determin outliers
  # previously we used 4x sd
  psm_df <- psm_df %>% mutate(isLow = ratioQC < meanRatio - nSD*sdRatio)
  psm_df <- psm_df %>% mutate(isHigh = ratioQC > meanRatio + nSD*sdRatio)
  psm_df <- psm_df %>% mutate(isOutlier = isLow | isHigh)
  # summary
  summary_df <- psm_df %>% group_by(bin) %>%
  	summarize(meanIntensity = mean(Intensity),
  		  meanRatio = mean(ratioQC,na.rm=TRUE),
  		  sdRatio = sd(ratioQC), 
  		  nOut = sum(isOutlier), .groups="drop")
  # FIXME: really this should be done for each experiment independently!
  outlier_psm <- psm_df %>% ungroup() %>% 
  	filter(isOutlier) %>% dplyr::select(PSM) %>% unlist() %>% unique()
  return(outlier_psm)
}

message("\nNumber of outlier PSM identified: ", length(outlier_psm))
