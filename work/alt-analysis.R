#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: Preprocessing and statistical analysis of Swip TMT proteomics.
#' authors: Tyler W A Bradshaw
#' ---

## INPUT data is in root/data/TMT_.zip.
# The zipped directory contains sample meta data and raw peptide data from PD:
zip_file = "TMT.zip" 
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

## OPTIONS:
fig_height = 5.0 # Default height of figures (in).
fig_width = 5.0 # Default width of figures (in).

FDR_alpha = 0.1 # FDR threshold for differential abundance.
BF_alpha = 0.05 # Bonferroni threshold for differential abundance.

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

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases 
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.

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
datadir <- file.path(rootdir, "data") # Key pieces of data saved as rda.
rdatdir <- file.path(rootdir, "rdata") # Temporary data files.
downdir <- file.path(rootdir, "downloads") # Misc downloads/temporary files.
tabsdir <- file.path(rootdir, "tables") # Output tables saved as excel files.

# Create project output directories if necessary.
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
entrez <- mgi_batch_query(uniprot)
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


#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

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
## Summarize to protein level.
#---------------------------------------------------------------------

message("\nSummarizing proteins as the sum of their peptides.")
proteins <- summarize_prot(filt_peptide)

# Perform SL normalization across all experiments (grouped by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")


#---------------------------------------------------------------------
## Perform IRS Normalization.
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
## Perform Sample Pool Normalization.
#---------------------------------------------------------------------

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
filt_protein <- filtProt(spn_protein,
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

# At this point there are no remaining missing values
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }


#---------------------------------------------------------------------
# combine the final normalized data and sample meta data
#---------------------------------------------------------------------

cols <- intersect(colnames(samples),colnames(filt_protein))
swip_tmt <- filt_protein %>% 
	left_join(samples,  by = cols) %>%
	filter(Treatment != "SPQC") %>%
	mutate(Genotype = Treatment) %>%
	mutate(Protein = Accession) %>% 
	mutate(Mixture = gsub("Exp","M", Experiment)) %>%
	mutate(BioFraction = Fraction) %>%
	mutate(Abundance = log2(Intensity)) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Intensity,Abundance)
irs_prot <- irs_protein %>% 
	left_join(samples,  by = cols) %>%
	filter(Treatment != "SPQC") %>%
	mutate(Genotype = Treatment) %>%
	mutate(Protein = Accession) %>% 
	mutate(Mixture = gsub("Exp","M", Experiment)) %>%
	mutate(BioFraction = Fraction) %>%
	mutate(Abundance = log2(Intensity)) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Intensity,Abundance)


proteins <- unique(swip_tmt$Protein)
results_list <- list()

library(doParallel)

doParallel::registerDoParallel(parallel::detectCores() - 1)

results_list <- foreach(prot = proteins) %dopar% {
  fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)
  df <- swip_tmt %>% 
	subset(Protein == prot) %>% 
	mutate(rel_Intensity = Intensity/sum(Intensity))
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, df, control = lmer_control)
  LT <- getContrast(fm,"Mutant","Control")
  result = lmerTestContrast(fm,LT) %>% 
  	     mutate(Contrast='Mutant-Control') %>% 
	     mutate(Protein = prot) %>%
	     unique()
  return(result)
}

df = dplyr::bind_rows(results_list)

df = df %>% mutate(FDR = p.adjust(Pvalue, method = "BH"))

sum(df$FDR<0.05)

data(washc_prots)
df %>% filter(Protein %in% washc_prots)

df$Symbol <- mapID(df$Protein,"uniprot","symbol")

fwrite(df, "foo.csv")


#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------
# Save key results.

# Save tidy_protein (final normalized protein in tidy format) as rda object. 
myfile <- file.path(datadir,"swip_tmt.rda")
save(swip_tmt,file=myfile,version=2)

myfile <- file.path(datadir,"irs_prot.rda")
save(irs_prot,file=myfile,version=2)

message("\nDone!")
