#!/usr/bin/env Rscript

# title: SwipProteomics
# description: SWIP TMT data preprocessing
# author: twab

## ---- INPUT 

# Load renv -- use renv::load NOT activate!
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE) 

# input data is in root/data/TMT.zip 
# TMT.zip contanis sample meta data and raw peptide data from PD
zip_file = "TMT.zip" 
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

FDR_alpha = 0.05 # FDR threshold for differential abundance
BF_alpha = 0.05 # Bonferroni threshold for differential abundance


## OUTPUT saved in root/.../tables:
# * Swip_TMT_Results.xlsx - an excel spreadsheet that contains:
#     - Sample meta data
#     - Gene/protein identifiers
#     - Raw peptide data
#     - Final normalized data (log2 transformed)
#     - Statistical results for intrafraction contrasts (F4-10)
#     - Statistical results for WT v MUT contrast

## OUTPUT data in root/data:
# Key datasets are saved in the data directory as read-only objects.
# * tmt_protein.rda

## OUTPUT for downstream analysis in root/rdata/:
# Input files/temporary files that are passed to other scripts including 
# the Python Leidenalg clustering script are saved in root/rdata.

## Order of data processing operations:
# * Load the raw peptide intensity data from PD 2.2
# * Sample loading normalization
# * Impute missing peptide values with KNN algorithm (k=10)
# * Examine QC peptide reproducibility - remove QC outliers
# * Summarize to protein level by summing all peptides for a protein
# * Protein-level sample loading normalization
# * IRS normalization - equalize protein measurements from different peptides
# * Sample pool normalization - assumes that mean of pooled QC technical 
#   replicates is equal to mean of all biological replicates
# * Protein level filtering - remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates)


## ---- prepare the renv

# Prepare the R workspace for the analysis

# Load required packages and functions
suppressPackageStartupMessages({
	library(dplyr)
	library(getPPIs)
	library(data.table)
	library(doParallel)
})

# Load project specific functions and data
suppressWarnings({ devtools::load_all(root) })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata") 
tabsdir <- file.path(root, "tables") 
downdir <- file.path(root, "downloads")

# Create project output directories if necessary
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }


## ---- Load the raw data and sample info

# Extract the raw TMT data from zipped file.
myfile <- file.path(datadir,zip_file)
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# Read TMT.csv data into R with data.table::fread
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_data)
peptides <- fread(myfile)

# Load sample information
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- fread(myfile)

# Format Cfg force column -- this is the Force in g's used to obtain 
# the subcellular fraction
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")


## ---- Map all Uniprot accession numbers to stable entrez IDs

# Map Uniprot IDs to Entrez using MGI batch query function
# Map Entrez IDs to gene symbols using getPPIs::getIDs()
message("\nCreating gene identifier map.")

# First, remove any non-mouse proteins from the data
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# Remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# collect all uniprot ids
uniprot <- unique(peptides$Accession)

# Map Uniprot IDs to Entrez using online MGI batch query function
# This takes a couple minutes because the function currently
# downloads the MGI data each time the function is called.
entrez <- mgi_batch_query(uniprot)
names(entrez) <- uniprot

# map any remaining missing ids by hand
message("Mapping missing IDs by hand.\n")
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(P05214=22144, P0CG14=214987, P84244=15078)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }

# Map entrez ids to gene symbols using twesleyb/getPPIs.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Create gene identifier mapping data.table
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")


## ---- Tidy the input data from Proteome Discover

# Convert PD df into tidy df
# Samples should contain the following columns:
# Treatment, Channel, Sample, Experiment
message("\nLoading raw data from Proteome Discover (PD2.2).")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,id.vars=cols)

# Annotate tidy data with additional meta data from samples
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


## ---- Perform sample loading normalization

# Perform sample normalization. Normalization is done for each 
# experiment independently (group by Experiment:Sample).
# NOTE: Grouping by Experiment, Channel does't work because 
# Sample TMT channels (e.g. 126N were used in different expirements).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))


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


#---------------------------------------------------------------------
## Summarize to protein level
#---------------------------------------------------------------------

message("\nSummarizing proteins as the sum of their peptides.")
proteins <- summarize_prot(filt_peptide)

# Perform SL normalization across all experiments (grouped by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")

# add additional meta data
sl_protein <- left_join(sl_protein,samples,by= intersect(colnames(sl_protein),colnames(samples)))


#---------------------------------------------------------------------

prot = swip
prot = sample(unique(irs_protein$Accession),1)
fx = log2(Intensity) ~ (1|Genotype) + (1|Fraction) + (1|Experiment)
lmer_args <- list()
lmer_args[["formula"]] <- fx
lmer_args[["data"]] <- sl_protein %>% filter(Treatment != "SPQC") %>% subset(Accession == prot)
fm <- do.call(lmerTest::lmer, lmer_args)
vp <- getVariance(fm)
pve <- vp/sum(vp)
round(pve,3)

#---------------------------------------------------------------------
## Perform IRS Normalization
#---------------------------------------------------------------------

# Equalize QC measurements between experiments. Adjusts protein 
# measurements of biological replicates simultaneously.
# Accounts for protein quantification by different peptides in 
# each experiment.
message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(sl_protein,controls="SPQC",robust=TRUE)

## Check, protein-wise QC measurements are now equal in a random protein.
## Generate list of QC proteins, and check the average of the three Exps.
#qc_proteins <- irs_protein %>% filter(Treatment == "SPQC") %>% 
#	group_by(Accession,Experiment) %>% 
#	dplyr::summarize(Treatment = unique(Treatment),
#		  "Mean(Intensity)" = mean(Intensity,na.rm=TRUE),
#		  .groups = "drop") %>% group_by(Accession) %>% group_split()
#message("\nIntra-experimental QC means are equal after IRS normalization:")
#knitr::kable(sample(qc_proteins,1))

#---------------------------------------------------------------------
## Perform Sample Pool Normalization
#---------------------------------------------------------------------

# QC samples were generated by pooling all biological replicates.
# Normalize mean of all QC samples to be equal to the mean of all 
# biological replicates.
message(paste("\nAccounting for experimental batch effect by",
	      "performing sample pool normalization."))
spn_protein <- normSP(irs_protein, pool=c("Control","Mutant"))


#---------------------------------------------------------------------

prot = sample(unique(spn_protein$Accession),1)
fx = log2(Intensity) ~ (1|Genotype) + (1|Fraction) + (1|Experiment)
lmer_args <- list()
lmer_args[["formula"]] <- fx
lmer_args[["data"]] <- sl_protein %>% filter(Treatment != "SPQC") %>% subset(Accession == prot)
fm <- do.call(lmerTest::lmer, lmer_args)
vp <- getVariance(fm)
pve <- vp/sum(vp)
round(pve,3)

#---------------------------------------------------------------------
## Protein level filtering
#---------------------------------------------------------------------

# Remove proteins that:
# * Were identified by a single peptide.
# * Contain too many (>50%) missing values.
# * Contain any missing QC values.
# * Have outlier protein measurements: 
#   mean(logRatio(replicates)) outside +/- nSD * mean(binned logRatio))
message(paste("\nFiltering proteins, this may take several minutes."))
filt_protein <- filtProt(spn_protein, controls="SPQC",nbins=5,nSD=4,summary=T)

# At this point there are no remaining missing values
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

#---------------------------------------------------------------------

## statistical testing with msstats
fx <- log2(Intensity) ~ 0 + Genotype:Fraction + (1|Experiment)
lmer_args <- list()
lmer_args[["data"]] <- filt_protein %>% subset(Accession == swip)
lmer_args[["formula"]] <- fx
fm <- lmerFit(lmer_args)
LT <- getContrast(fm,"MUT","WT")
lmerTestContrast(fm,LT)

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# fit protein-wise models
proteins <- unique(filt_protein$Accession)
fit_list <- foreach (protein = proteins) %dopar% {
	lmer_args <- list()
	lmer_args[["formula"]] <- fx
	lmer_args[["data"]] <- filt_protein %>% subset(Accession == protein)
	fm <- lmerFit(lmer_args)
	return(fm)
} #EOL
names(fit_list) <- proteins


# assess protein-wise contrasts with lmerTestContrast
results_list <- foreach (protein = proteins) %dopar% {j
	fm <- fit_list[[protein]]
	lmerTestContrast(fm,LT) %>% 
		mutate(Contrast="Mutant-Control") %>% unique()
} #EOL

names(results_list)

results_df <- bind_rows(results_list,.id="Protein")

#---------------------------------------------------------------------

tmt_prot <- filt_protein %>% filter(Treatment != "SPQC")
colnames(tmt_prot)[colnames(tmt_prot) == "Experiment"] <- "Mixture"
colnames(tmt_prot)[colnames(tmt_prot) == "Treatment"] <- "Genotype"
colnames(tmt_prot)[colnames(tmt_prot) == "Accession"] <- "Protein"
tmt_prot <- tmt_prot %>% mutate(Abundance = log2(Intensity)) %>%
	mutate(Genotype = samples$Genotype[match(Sample,samples$Sample)]) %>%
	mutate(BioFraction = samples$Fraction[match(Sample,samples$Sample)]) %>%
	mutate(Condition = interaction(Genotype,BioFraction))
tmt_prot$Intensity <- NULL
tmt_prot$Sample <- NULL

# cast the data into a matrix
dm <- tmt_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Condition, value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

# check for na
idx <- apply(dm,1,function(x) any(is.na(x)))
stopifnot(sum(idx)==0)

# calc corr mat
adjm <- cor(t(dm),method="pearson")

# network enhancement
ne_adjm <- neten::neten(adjm)
