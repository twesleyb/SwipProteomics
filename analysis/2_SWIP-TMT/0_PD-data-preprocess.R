#!/usr/bin/env Rscript

# title: SwipProteomics
# description: prepare the PD PSM level data for analysis by MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## Input ----------------------------------------------------------------------
# specify the projects root directory:
root <- "~/projects/SwipProteomics"

# the PSM data is in root/rdata:
# NOTE: the raw data is 75 mb and is too large to be managed by git.
# After preprocessing, the PSM-level data is saved as pd_psm.rda in
# root/data.
input_dir <- "rdata/PSM.zip"

# PSM.zip contains:
input_psm <- "PSM/5359_PSM_Report.xlsx" # the data exported from PD
input_samples <- "PSM/5359_Sample_Report.xlsx" # MS run and sample info


## Options --------------------------------------------------------------------
save_rda <- TRUE # save key R objects?


## Output ---------------------------------------------------------------------
# saved in root/data

# * gene_map - data.table mapping UniProt Accession to Entrez IDs and gene
#       Symbols.
# * pd_psm - PSM data from Proteome Discoverer reformatted for MSstatsTMT.
#       This file is just less than 50 mb.
# * pd_annotation - the annotation object for PDtoMSstatsTMTformat().
#       This maps individual MS runs to Spectrum.File and other sample
#       information
# * msstats_contrasts - a matrix specifying all pairwise intrafraction contrasts
#       between Control and SWIP P1019R homozygous Mutant mice.


## Functions -----------------------------------------------------------------
# functions utilized herein

squote <- function(string) {
  # wrap a string in singular quotation marks (')
  return(paste0("'", string, "'"))
}


mkdir <- function(..., warn = TRUE, report = FALSE) {
  # create a new directory
  newdir <- file.path(...)
  if (warn & dir.exists(newdir)) {
    warning("dir exists")
  } else if (!dir.exists(newdir)) {
    dir.create(newdir)
    if (report) {
      message("\nCreated ", newdir, ".")
    }
  }
}


fix_colname <- function(df, old_colname, new_colname) {
  # change a column's name in a data.frame
  colnames(df)[which(colnames(df) == old_colname)] <- new_colname
  return(df)
}


munge1 <- function(x) { # x = samples$ConditionFraction
  # function to reformat the data
  # extract sample 'Condition' annotation from ConditionFraction
  paste0("F", as.numeric(sapply(strsplit(x, "Control|Mutant|SPQC"), "[", 2)))
}


munge2 <- function(x) { # x = samples$ConditionFraction
  # function to reformat the data
  # extract sample 'Fraction' annotation from ConditionFraction
  sapply(strsplit(x, "[0-9]{1,2}"), "[", 1)
}


munge3 <- function(x) { # x = samples$Experiment
  # coerce Experiment to replicate ID = R#
  paste0("R", as.numeric(as.factor(samples$Experiment)))
}


reformat_cols <- function(raw_pd) {
  # make columns look like MSstats by replacing special characters with '.'
  # replace special characters in column names with "."
  chars <- c(" ", "\\[", "\\]", "\\:", "\\(", "\\)", "\\/", "\\+", "\\#", "\\-")
  new_names <- gsub(paste(chars, collapse = "|"), ".", colnames(raw_pd))
  colnames(raw_pd) <- new_names
  # add 'X' if starts column name starts with ".."
  colnames(raw_pd) <- gsub("^\\.\\.", "X..", colnames(raw_pd))
  # return the reformatted data
  return(raw_pd)
}


## Prepare the working environment --------------------------------------------

# project directories
datadir <- file.path(root, "data")
mkdir(datadir, warn = FALSE)
rdatdir <- file.path(root, "rdata")
mkdir(rdatdir, warn = FALSE)
downdir <- file.path(root, "downloads")
mkdir(downdir, warn = FALSE)


# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root, quiet = TRUE)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  suppressWarnings({
    library(getPPIs)
  }) # FIXME: annoying warnings!
  library(data.table)
  library(MSstatsTMT) # twesleyb/MSstatsTMT
})

# load functions in root/R
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})


## unzip the data in root/data ------------------------------------------------

# unzip data into downloads
unzip(file.path(root, input_dir), exdir = downdir)
message("\nUnzipped ", squote(input_dir), " into ", squote(downdir), ".")


## read PSM data from excel ---------------------------------------------------

# read PSM-level data exported from PD as an excel worksheet
# NOTE: this takes several minutes!

# load the data
myfile <- file.path(downdir, input_psm)
raw_pd <- readxl::read_excel(myfile, progress = FALSE)

# re-format PSM data for MSstatsTMT
raw_pd <- reformat_cols(raw_pd) # changes colnames to match what MSstats expects


## load sample data -----------------------------------------------------------
# recieved this excel spreadsheet from GW--I think he exported from PD

# pass meaningful colnames to read_excel
col_names <- c(
  "Sample", "Mixture", "MS.Channel", "drop",
  "Channel", "Proteomics ID", "ConditionFraction", "Experiment"
)
myfile <- file.path(downdir, input_samples)
samples <- readxl::read_excel(myfile, col_names = col_names)

# clean-up
mydir <- file.path(downdir, tools::file_path_sans_ext(basename(input_dir)))
unlink(mydir, recursive = TRUE)
message("\nRemoved ", squote(mydir), ".")


## re-format sample metadata annotations for MSstats ---------------------------

# remove un-needed col
samples$drop <- NULL

# Munge 'ConditionFraction' column to 'BioFraction' and 'BioCondition'
# NOTE: BioFraction is the subcellular fraction not MSstats 'Fraction'
# NOTE: BioCondition is the treatment 'Condition'
samples$BioFraction <- munge1(samples$ConditionFraction)
samples$Genotype <- munge2(samples$ConditionFraction)

# remove un-needed col
samples$ConditionFraction <- NULL

# this is how MSstatsTMT needs 'Condition' for intra-fraction contrasts:
# >>> BioCondition.BioFraction e.g. Control.F10
condition <- as.character(interaction(samples$Genotype, samples$BioFraction))
condition[grepl("SPQC", condition)] <- "Norm"
samples$Condition <- condition

# clean-up 'Mixture' column by replacing F with M
samples$Mixture <- gsub("F", "M", samples$Mixture)

# FIXME: how should BioReplicate be defined?
# for intrafraction contrasts, 'BioReplicate' should be ...
# add R# to indicate the replicate %eg% Control.F4.R2
biorep <- as.character(interaction(samples$Experiment, samples$Condition))
biorep[grep("Norm", biorep)] <- "Norm"
samples$BioReplicate <- biorep

# Subject
samples$Subject <- samples$BioReplicate


## drop contaminant proteins ---------------------------------------------------

# drop PSM that are mapped to multiple proteins
idx_drop <- grepl(";", raw_pd$Master.Protein.Accessions)
filt_pd <- raw_pd[!idx_drop, ]

# keep mouse proteins and drop ig prots
idx_keep <- grepl("OS=Mus musculus", filt_pd$Protein.Descriptions)
ig_prots <- c(
  "P01631", "P01646", "P01665", "P01680", "P01746", "P01750",
  "P01786", "P01864", "P01878", "P03975", "P06330", "P03987"
)

# drop human keratins, proteins without genes, pig trypsin, a predicted gene
misc_drop <- c(
  "P13647", "P13645", "P08779", "P04264", "P02533", "Q04695",
  "Q7Z794", "P01626", "P01637", "P00761", "P10404", "P02538", "Q3ZRW6",
  "P04940", "Q6KB66", "P01723", "Q80WG5"
)

# filter data
idx_drop1 <- filt_pd$Master.Protein.Accessions %in% ig_prots
idx_drop2 <- filt_pd$Master.Protein.Accessions %in% misc_drop
filt_pd <- filt_pd[idx_keep & !idx_drop1 & !idx_drop2, ]

# collect all uniprot accession ids
uniprot <- unique(filt_pd$Master.Protein.Accessions)


## create gene map by mapping uniprot to entrez -------------------------------
# NOTE: this may take several minutes

# all entrez gene ids
entrez <- suppressWarnings({
  mgi_batch_query(uniprot, quiet = TRUE)
})
names(entrez) <- uniprot

# Map any remaining missing IDs by hand.
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(
  P05214 = 22144,
  P0CG14 = 214987,
  P84244 = 15078,
  P05214 = 22144
) # tub3a
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
if (!check) {
  stop("Unable to map all UniprotIDs to Entrez.")
}

# Map entrez ids to gene symbols using twesleyb/getPPIs.
symbols <- getPPIs::getIDs(entrez,
  from = "entrez",
  to = "symbol", species = "mouse"
)

# check there should be no missing gene symbols
if (any(is.na(symbols))) {
  stop("Unable to map all Entrez to gene Symbols.")
}
# Create gene identifier mapping data.table.
gene_map <- data.table(
  uniprot = names(entrez),
  entrez = entrez,
  symbol = symbols
)
gene_map$id <- paste(gene_map$symbol, gene_map$uniprot, sep = "|")


## map Spectrum.Files to MS.Channels ------------------------------------------
# Basically, its just complicated. There were 3x 16-plex TMT experiments.
# In each, the concatenated TMT mixture was fractionated into 12 fractions in
# order to increase analytical depth. Therefore, there were 12 * 3 = 36 mass
# spectrometry runs. Each Run is recorded as a 'Spectrum File' by Proteome
# Discover. Note that each MS run cooresponds to the measurment of all 16
# samples and therefore the total number of TMT channels for which we
# have made measurements is 16 x 3 Experiments = 48. In other words, a single
# 'Spectrum.File' cooresponds to 12x MS.Runs and 16x Samples.

all_files <- filt_pd$Spectrum.File

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_" to extract experiment identifiers
exp_files <- lapply(
  split(all_files, sapply(strsplit(all_files, "_"), "[", 1)),
  unique
)
files_dt <- data.table(
  "Experiment" = rep(names(exp_files),
    times = sapply(exp_files, length)
  ),
  "Run" = unlist(exp_files)
)

# add Fraction annotation
files_dt$Fraction <- unlist({
  sapply(exp_files, function(x) as.numeric(as.factor(x)), simplify = FALSE)
})

# collect all MS.Channels, grouped by Experiment
all_channels <- samples$MS.Channel
exp_channels <- split(all_channels, sapply(strsplit(all_channels, "_"), "[", 1))

# as a dt
exp_dt <- data.table(
  "Experiment" = rep(names(exp_channels),
    times = sapply(exp_channels, length)
  ),
  "MS.Channel" = unlist(exp_channels)
)


## build annotation file for MSstats -------------------------------------------
# the annotation data.frame passed to MSstats requires the following columns:
# * Run - indicates which channel within a Spectrum.File; this should match
#     Spectrum.File in raw_pd.
# * Fraction - a TMT mixture may be fractionated into multiple fractions to
#     increase analytical depth; column 'Fraction' indicates the MS fraction,
#     it may be all 1.
# * TechRepMixture - was a TMT mixture analyzed in technical replicate?
#     all 1 indicates no replicates of a mixture.
# * Mixture - concatenation of multiple TMT labeled samples - an MS experiment.
# * Channel - the TMT labels/channels used e.g. 126N, 134N
# * BioReplicate - indicates an individual Biological subject
# * Condition - indicates the treatment condition e.g. WT, MUT or SPQC (Norm)
# NOTE: MSstats expects the Norm(alization) condition to be 'Norm'.

# create annotation_dt from Spectrum.Files and MS.Runs
# add additional freatures from samples
annotation_dt <- left_join(files_dt, exp_dt, by = "Experiment")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$BioFraction <- samples$BioFraction[idx]
annotation_dt$TechRepMixture <- rep(1, length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
annotation_dt$Channel <- samples$Channel[idx]
annotation_dt$BioReplicate <- samples$BioReplicate[idx]

# Remove un-needed cols
annotation_dt$"MS.Channel" <- NULL


## build contrast_matrix ----------------------------------------------------
# create MSstatsTMT contrasts for intrafraction comparisons

# define all intrafraction comparisons:
comp <- paste(paste("Mutant", paste0("F", seq(4, 10)), sep = "."),
  paste("Control", paste0("F", seq(4, 10)), sep = "."),
  sep = "-"
)

# create a contrast matrix for given comparisons
conditions <- unique(annotation_dt$Condition)
conditions <- conditions[conditions != "Norm"] 

# utilizes internal function made available by my fork to generate a contrast
# matrix for all pairwise comparisions defined by comp
all_contrasts <- MSstatsTMT::makeContrast(groups = conditions)

# subset contrast matrix, keep pairwise contrasts of interest
biof <- sapply(strsplit(rownames(all_contrasts), "\\.|-"), "[",
  c(2, 4),
  simplify = FALSE
)
idx <- sapply(biof, function(x) x[1] == x[2])

# subset and flip sign of the comparisons: Mutant - Control
msstats_contrasts <- -1 * all_contrasts[idx, ]

# flip rownames
namen <- rownames(msstats_contrasts)
new_names <- paste(sapply(strsplit(namen,"-"),"[",2),
		   sapply(strsplit(namen,"-"),"[",1),sep="-")
rownames(msstats_contrasts) <- new_names


## Create 'Mutant-Control contrast ---------------------------------------------

# add a contrast vector specifying the mutant-control comparison to
# msstats_contrasts matrix (see 0_*.R)

mut_vs_control <- matrix(c(
  -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7, -1 / 7,
  1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7, 1 / 7
), nrow = 1)
row.names(mut_vs_control) <- "Mutant-Control"
colnames(mut_vs_control) <- colnames(msstats_contrasts)


## Save outputs to file ---------------------------------------------------------

# save WASHC4's uniprot ID
swip <- gene_map$uniprot[gene_map$symbol == "Washc4"]
myfile <- file.path(root, "data", "swip.rda")
save(swip, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)

# gene_map - gene identifiers for all proteins
myfile <- file.path(datadir, "gene_map.rda")
save(gene_map, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)

# pd_annotation - input annotations for MSstatsTMT
pd_annotation <- annotation_dt
myfile <- file.path(datadir, "pd_annotation.rda")
save(pd_annotation, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)

# pd_psm - the raw PSM data--input for MSstatsTMT
pd_psm <- filt_pd
myfile <- file.path(datadir, "pd_psm.rda")
save(pd_psm, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)

# msstats_contrasts - all pairwise intraBioFraction comparisons
myfile <- file.path(datadir, "msstats_contrasts.rda")
save(msstats_contrasts, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)

# mut_vs_control - the 'Mutant-Control' comparison
myfile <- file.path(datadir, "mut_vs_control.rda")
save(mut_vs_control, file = myfile, version = 2)
message(
"\nSaved ", squote(basename(myfile)), " in ",
squote(dirname(myfile)), "."
)
