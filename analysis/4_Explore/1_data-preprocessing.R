#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: Preprocessing and statistical analysis of Swip TMT proteomics.
# author: Tyler W A Bradshaw

## INPUT data in root/data/TMT_.zip -------------------------------------------
# The zipped directory contains sample meta data and raw peptide data from PD:
zip_file = "TMT.zip" 
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

## OPTIONS ---------------------------------------------------------------------
fig_height = 5.0 # Default height of figures (in).
fig_width = 5.0 # Default width of figures (in).

oldham_threshold = 2.5 # Sample connectivity threshold for detecting outliers.
# NOTE: No sample outliers were detected.

FDR_alpha = 0.1 # FDR threshold for differential abundance.
BF_alpha = 0.05 # Bonferroni threshold for differential abundance.
logFC_threshold = c(lwr=log2(0.8), upr=log2(1.2)) # FC threshold.

set.seed = as.numeric(Sys.time()) # seed for random operations.

## OUTPUT saved in root/.../tables --------------------------------------------
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

## DESCRIPTION ----------------------------------------------------------------

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

## Functions ------------------------------------------------------------------
## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases 
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly. In addition to 
## any functions declared in this chunk functions from root/R are utilitzed.

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

## Prepare the R environment --------------------------------------------------
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
datadir <- file.path(rootdir, "data") # Key pieces of data saved as rda.
rdatdir <- file.path(rootdir, "rdata") # Temporary data files.
downdir <- file.path(rootdir, "downloads") # Misc downloads/temporary files.
tabsdir <- file.path(rootdir,"manuscript","tables") # Output tables saved as excel files.

# Create output directories if necessary.
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

## Load the raw data and sample info ------------------------------------------

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

## Map all Uniprot accession numbers to stable entrez IDs ---------------------

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

## Tidy-up the input data from Proteome Discover ------------------------------

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

## Perform sample loading normalization ---------------------------------------

# Perform sample normalization. Normalization is done for each 
# experiment independently (group by Experiment:Sample).
# NOTE: Grouping by Experiment, Channel does't work because 
# Sample TMT channels (e.g. 126N were used in different expirements).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))

# Check the data before SL:
message("\nPeptide data before sample loading normalization:")
check_SL(tidy_peptide)

# Check the data after SL:
message("\nPeptide data after sample loading normalization:")
check_SL(sl_peptide) # Equal within an experiment.

## Impute missing peptide values -----------------------------------------------

# Impute missing peptide values with k-nearest neighbors (KNN) algorithm.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
# NOTE: KNN imputing is done for each experimental group seperately.
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC",quiet=FALSE) 

## Examine reproducibility of QC measurements ----------------------------------
# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated. 
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, measuremnts that are outside 
# +/- 4x standard deviations from the bin's mean are removed. 

message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide,controls="SPQC",quiet=FALSE)

## Summarize to protein level -------------------------------------------------

# Protein intensity is just the sum of all peptides for that protein.
message("\nSummarizing proteins as the sum of their peptides.")
proteins <- summarize_prot(filt_peptide)

# Perform SL normalization across all experiments (grouped by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")

# Column sums are equal across all experiments:
message("\nProtein data after sample loading normalization:")
check_SL(sl_protein)

## Perform IRS Normalization --------------------------------------------------

# Equalize QC measurements between experiments. Adjusts protein 
# measurements of biological replicates simultaneously.
# Accounts for protein quantification by different peptides in 
# each experiment.
message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(sl_protein,controls="SPQC",robust=TRUE)

# Check, protein-wise QC measurements are now equal in a random protein.
# Generate list of QC proteins, and check the average of the three Exps.
qc_proteins <- irs_protein %>% filter(Treatment == "SPQC") %>% 
	group_by(Accession,Experiment) %>% 
	dplyr::summarize(Treatment = unique(Treatment),
		  "Mean(Intensity)" = mean(Intensity,na.rm=TRUE),
		  .groups = "drop") %>% group_by(Accession) %>% group_split()
message("\nIntra-experimental QC means are equal after IRS normalization:")
knitr::kable(sample(qc_proteins,1))

## Perform Sample Pool Normalization ------------------------------------------

# QC samples were generated by pooling all biological replicates.
# Normalize mean of all QC samples to be equal to the mean of all 
# biological replicates.
message(paste("\nAccounting for experimental batch effect by",
	      "performing sample pool normalization."))
spn_protein <- normSP(irs_protein,pool=c("Control","Mutant"))

# Check, the mean of sample pool and biological replicates are now equal.
message(paste("\nThe mean of biological replicates and",
	      "pooled QC samples are now equal:"))
protein_list <- spn_protein %>% group_by(Accession)  %>% group_split()

# Drop any proteins with NA.
idx <- sapply(protein_list,function(x) any(is.na(x)))
protein_list <- protein_list[!idx]

# Get a random protein's data as an example.
protein <- sample(protein_list,1)[[1]]
protein$Sample_Pool <- protein$Treatment == "SPQC"
df <- protein %>% group_by(Sample_Pool) %>% 
	dplyr::summarize(Accession = unique(Accession),
		  Treatement=paste(unique(Treatment),collapse=" + "),
		  "Mean(Intensity)"=mean(Intensity), .groups="drop")
knitr::kable(df)

rm(list=c("df","protein_list","idx","protein"))

## Protein level filtering ----------------------------------------------------
# Remove proteins that:
# * Were identified by a single peptide.
# * Contain too many (>50%) missing values.
# * Contain any missing QC values.
# * Have outlier protein measurements: 
#   mean(logRatio(replicates)) outside +/- nSD * mean(binned logRatio))
message(paste("\nFiltering proteins, this may take several minutes."))
filt_protein <- filtProt(spn_protein,
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

# At this point there are no remaining missing values.
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

rm(list="check")

## Protein level imputing -----------------------------------------------------

# If necessary, protein level imputing can be performed.
# Proteins with more than 50% missing values are ignored.

message(paste("\nImputing missing protein values."))
imputed_protein <- imputeKNNprot(filt_protein,ignore="SPQC")

## Identify Sample outliers --------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(filt_protein)
outlier_samples <- c(names(zK)[zK < -1 * oldham_threshold],
		     names(zK)[zK > +1 * oldham_threshold])

# There are no sample outliers.
check <- length(outlier_samples) == 0
if (check) { 
	message("\nNo sample level outliers were identified.")
} else { 
	stop("Why are there outlier samples?") 
}



## Create Summarized experiment object with colData and meta data
colData <- as.data.frame(cbind(sample,group))
rownames(colData) <- sample
colData

# Convert to data matrix
dm <- as.matrix(df)

# Insure 0 values are NA
dm[dm==0] <- NA

# Create SE object
library(SummarizedExperiment)
data_se <- SummarizedExperiment(assays=list(tmt=dm), colData=colData)

# Add meta data
metadata <- list()
metadata$sample <- "sample"
metadata$group <- "group"
metadata(data_se) <- metadata
data_se

# Create Normalyzer DE object
library(NormalyzerDE)

jobname <- "Combined_Commmon_prots"
jobDir <- getwd()

normObj <- getVerifiedNormalyzerObject(jobname,data_se)

# Calc normalizations
normResults <- normMethods(normObj)

###############################################
## SL Normalization 

# Extract log2 data
log2.dm<-normResults@normalizations[["log2"]]
head(log2.dm)

# Unlog
dm <- 2.^log2.dm
head(dm)

# Column sums
col_sums <- colSums(dm,na.rm=TRUE)
head(col_sums)

# Avg of column sums
avg <- mean(col_sums)
avg

# Scaling factor
factors <- avg/col_sums
factors

# Normalyze by scaling factor
sl.dm <- sweep(dm, 2, factors, FUN = "*")

# Check column sums should be equal
colSums(sl.dm,na.rm=TRUE)


# Log transform
sl.dm <- log2(sl.dm)

## add to normResults object:
normResults@normalizations$SL<- sl.dm

###############################################
## SL+IRS Normalization
qc_cols <- grep("QC",colnames(sl.dm))
qc_dm <- sl.dm[,qc_cols]
av1 <- rowMeans(qc_dm[,1:3],dims=1)
av2 <- rowMeans(qc_dm[,4:6],dims=1)
av3 <- rowMeans(qc_dm[,7:9],dims=1)
av4 <- rowMeans(qc_dm[,10:12],dims=1)
QC_avg <- as.matrix(cbind(av1,av2,av3))
avg_QC_avg <- rowMeans(QC_avg,dims=1) 

# Calculate factor as group average / global average 
f1 <- as.matrix(av1/avg_QC_avg)
f2 <- as.matrix(av2/avg_QC_avg)
f3 <- as.matrix(av3/avg_QC_avg)
f4 <- as.matrix(av4/avg_QC_avg)

# Matrix of factors
F1 <- matrix(rep(f1,11), ncol = 11)
F2 <- matrix(rep(f2,11), ncol = 11)
F3 <- matrix(rep(f3,11), ncol = 11)
F4 <- matrix(rep(f4,11), ncol = 11)
factors.dm <- cbind(F1,F2,F3,F4)

# Normalize cleandf.Imputed by calculated factors
sl.irs.dm <- sl.dm/factors.dm

# Check, average of QCs should be equal
qc_dm <- sl.irs.dm[,qc_cols]
av1 <- rowMeans(qc_dm[,1:3],dims=1)
av2 <- rowMeans(qc_dm[,4:6],dims=1)
av3 <- rowMeans(qc_dm[,7:9],dims=1)
av4 <- rowMeans(qc_dm[,10:12],dims=1)
QC_avg <- as.matrix(cbind(av1,av2,av3))
head(QC_avg)
# GOOD 

## add to normResults object:
normResults@normalizations$SL.IRS <- sl.irs.dm

###############################################
## VSN+IRS Normalization 

# Or Extract VSN normalized data
vsn.dm<-normResults@normalizations[["VSN"]]
head(vsn.dm)

qc_cols <- grep("QC",colnames(vsn.dm))
qc_dm <- vsn.dm[,qc_cols]
av1 <- rowMeans(qc_dm[,1:3],dims=1)
av2 <- rowMeans(qc_dm[,4:6],dims=1)
av3 <- rowMeans(qc_dm[,7:9],dims=1)
av4 <- rowMeans(qc_dm[,10:12],dims=1)
QC_avg <- as.matrix(cbind(av1,av2,av3))
avg_QC_avg <- rowMeans(QC_avg,dims=1) 

# Calculate factor as group average / global average 
f1 <- as.matrix(av1/avg_QC_avg)
f2 <- as.matrix(av2/avg_QC_avg)
f3 <- as.matrix(av3/avg_QC_avg)
f4 <- as.matrix(av4/avg_QC_avg)

# Matrix of factors
F1 <- matrix(rep(f1,11), ncol = 11)
F2 <- matrix(rep(f2,11), ncol = 11)
F3 <- matrix(rep(f3,11), ncol = 11)
F4 <- matrix(rep(f4,11), ncol = 11)
factors.dm <- cbind(F1,F2,F3,F4)

# Normalize cleandf.Imputed by calculated factors
normdf <- vsn.dm/factors.dm

## add to normResults object:
normResults@normalizations$VSN.IRS<- normdf

###############################################
## Recalculate performance metrics 
normResultsWithEval <- analyzeNormalizations(normResults)

# Generate plots.
generatePlots(normResultsWithEval, jobDir)

# Calculate performance metrics
normResultsWithEval <- analyzeNormalizations(normResults)

# Generate PDF output of evaluation plots. 
generatePlots(normResultsWithEval, jobDir)

#################################################################
## Evaluate differential expression using limma 

## Reformat for DEP 
#  Add ID and names columns
mergedf <- as.data.frame(normResults@normalizations$SL.IRS)
rawdf <- as.data.frame(normResults@normalizations$log2)



#################################################################################################
## Create SE object
library("DEP")
DEP_format <- function(data_in){
  data_in$ID <- rownames(data_in)
  data_in$name <- rownames(data_in)
  head(mergedf)
 TMT_columns <- 1:(ncol(data_in)-2)
 col.names <- colnames(data_in)[TMT_columns]
 col.names <- gsub("HET","HET.",col.names)
 col.names <- gsub("WT","WT.",col.names)
 col.names <- gsub("KO","KO.",col.names)
 col.names <- gsub("QC","QC.",col.names)
 col.names <- gsub("\\..*","",col.names)
 n <- (dim(data_in)[2]-2)/11
 replicate <- rep(c(1,2,3,1,2,3,4,1,2,3,4),n)
 col.names <- paste(col.names,replicate,sep=".")
 colnames(data_in)[TMT_columns] <- col.names
 return(data_in)
}

mergedf <- DEP_format(mergedf)
rawdf <- DEP_format(rawdf)

# Create SE object
data_raw <- make_se_parse(rawdf, TMT_columns)
data_se <- make_se_parse(mergedf, TMT_columns)

# Plot normalization
plot_normalization(data_raw,data_se)

# Filter rows with more than 50% missing values
data_filt <- filter_proteins(data_se, "fraction", min = 0.5)

# Examine remaining missing values:
plot_detect(data_filt)

# The distribution of missing values is shifted to the right. 
# MAR random values should be imputed with "knn" or "MLE"
# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp <- impute(data_filt, fun = "knn", rowmax = 0.5)
#data_imp <- impute(data_filt, fun = "MLE")

# Differential expression testing...
contrasts <- c("Syngap1HET._vs_Syngap1WT.","Shank3KO._vs_Shank3WT.", "Shank2KO._vs_Shank2WT.","Ube3aKO._vs_Ube3aWT.")
data_diff <- test_diff(data_imp, type = "manual", test = contrasts)

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

## Map uniprot IDs to gene names
# Get Uniprot_IDs
Uniprot_IDs <- data_results$ID
head(Uniprot_IDs)

# Map Uniprot IDs to Gene names 
library(AnnotationDbi)
library(org.Mm.eg.db)
symbol <- mapIds(org.Mm.eg.db,keys=Uniprot_IDs,column="SYMBOL", keytype="UNIPROT", multiVals="first")
data_results$name <- symbol

# write to table
write.csv(data_results,"Cortex_TMT_DEP_all_results.csv")
