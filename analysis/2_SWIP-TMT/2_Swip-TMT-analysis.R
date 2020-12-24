#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw

## ---- input

# project's root directory
root = "~/projects/SwipProteomics"

# INPUT data is in root/data/
zip_file = "TMT.zip"
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

# FDR threshold for differential abundance
FDR_alpha = 0.05


# ---- functions

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.


# ---- Prepare the workspace

# Prepare the R workspace for the analysis.

# Load renv -- use renv::load NOT activate!
renv::load(root, quiet=TRUE)

# Load required packages and functions
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(getPPIs) # twesleyb/getPPIs for mouse PPIs
	library(geneLists) # twesleyb/geneLists for gene mapping fun
	library(data.table) # for working with tables
	library(doParallel) # for parallel processing
})

# load project specific functions and data
devtools::load_all(root, quiet=TRUE)

# project directories:
datadir <- file.path(root, "data") # Key pieces of data saved as rda
rdatdir <- file.path(root, "rdata") # Temporary data files
tabsdir <- file.path(root, "tables") # Output tables saved as excel files
downdir <- file.path(root, "downloads") # Misc downloads/temporary files

# create project output directories if necessary
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }


## ---- Load the raw data and sample info

# Extract the raw TMT data from zipped file
myfile <- file.path(datadir, zip_file)
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# Read TMT.csv data into R with data.table::fread
myfile <- file.path(downdir, tools::file_path_sans_ext(zip_file), input_data)
peptides <- fread(myfile)

# Load sample information
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- fread(myfile)

# format cfg force column -- this is the Force in g's used to obtain
# the subcellular fraction.
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")


## ---- map all Uniprot accession numbers to stable entrez IDs

message("\nCreating gene identifier map.")

# first, remove any non-mouse proteins from the data
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# collect all Uniprot IDs
uniprot <- unique(peptides$Accession)

# map Uniprot IDs to Entrez using online MGI batch query function
entrez <- geneLists::queryMGI(uniprot)
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
# NOTE: getIDs is just an easier to use wrapper around AnnotationDbi mapIDs
# NOTE: you need the org.Mm.eg.db package for mapping mouse genes
gene_symbols <- geneLists::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Create gene identifier mapping data.table.
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = gene_symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")


## ---- Tidy-up the input data from Proteome Discover

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


##
sl_prot <- sl_peptide %>%
	group_by(Accession, Sample) %>%
	summarize(Intensity = sum(Intensity,na.rm=TRUE),.groups="drop") %>%
	filter(!is.na(Intensity)) %>% filter(Intensity != 0) %>%
	mutate(Abundance = log2(Intensity)) %>%
	left_join(samples, by="Sample") %>%
	filter(Treatment != "SPQC") %>%
	mutate(Genotype = Treatment) %>%
	mutate(Protein = Accession) %>%
	mutate(Mixture = gsub("Exp","M", Experiment)) %>%
	mutate(BioFraction = Fraction) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
	mutate(Subject = as.numeric(interaction(Mixture,Genotype))) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Subject,Abundance)
myfile = file.path(root,"rdata","sl_prot.rda")
save(sl_prot,file=myfile,version=2)


## ---- Perform sample loading normalization

# Perform sample normalization. Normalization is done for each
# experiment independently (group by Experiment:Sample).
# NOTE: Grouping by Experiment, Channel does't work because
# Sample TMT channels (e.g. 126N were used in different expirements).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))


## SAVE
cols <- intersect(colnames(samples),colnames(sl_peptide))
sl_prot <- sl_peptide %>%
	left_join(samples,  by = cols) %>%
	filter(Treatment != "SPQC") %>%
	mutate(Genotype = Treatment) %>%
	mutate(Protein = Accession) %>%
	mutate(Mixture = gsub("Exp","M", Experiment)) %>%
	mutate(BioFraction = Fraction) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
	mutate(Subject = as.numeric(interaction(Mixture,Genotype))) %>%
	group_by(Protein, Mixture, Genotype, BioFraction) %>%
	mutate(Abundance = log2(sum(Intensity,na.rm=TRUE))) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Subject,Abundance)
myfile = file.path(root,"rdata","sl_prot.rda")
save(sl_prot,file=myfile,version=2)


quit()

## ---- Impute missing peptide values

# Impute missing peptide values with k-nearest neighbors (KNN) algorithm.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
# NOTE: KNN imputing is done for each experimental group seperately.
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC", quiet=FALSE)


## ---- Examine reproducibility of QC measurements

# Assess reproducibility of QC measurements and remove QC peptide measurments
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of SPQC tech replicates is calculated.
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, ratios that are outside
# +/- 4x standard deviation from the bin's mean (should be centered at 0) are removed

message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide,controls="SPQC", quiet=FALSE)


## ---- Summarize to protein level

message("\nSummarizing proteins as the sum of their peptides.")
proteins <- sumProt(filt_peptide)

# Perform SL normalization across all experiments (grouped by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")


## ---- Perform IRS Normalization

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


## ---- Perform Sample Pool Normalization.

# QC samples were generated by pooling all biological replicates.
# Normalize mean of all QC samples to be equal to the mean of all
# biological replicates.
message(paste("\nAccounting for experimental batch effect by",
	      "performing sample pool normalization."))
spn_protein <- normSP(irs_protein,pool=c("Control","Mutant"))

# Check, the mean of sample pool and biological replicates are now equal.
message(paste("\nThe mean of biological replicates and",
	      "pooled QC samples are now equal:"))
protein_list <- spn_protein %>% group_by(Accession) %>% group_split()

# Drop any proteins with NA.
idx <- sapply(protein_list,function(x) any(is.na(x)))
protein_list <- protein_list[!idx]

# Get a random protein's data as an example
protein <- sample(protein_list,1)[[1]]
protein$Sample_Pool <- protein$Treatment == "SPQC"
df <- protein %>% group_by(Sample_Pool) %>%
	dplyr::summarize(Accession = unique(Accession),
		  Treatement=paste(unique(Treatment),collapse=" + "),
		  "Mean(Intensity)"=mean(Intensity), .groups="drop")
knitr::kable(df)


## ---- Protein level filtering

# Remove proteins that:
# * Were identified by a single peptide.
# * Contain too many (>50%) missing values.
# * Contain any missing QC values.
# * Have outlier protein measurements:
#       mean(logRatio(replicates)) outside +/- nSD * mean(binned logRatio))

#FIXME: filtProt is slow!

message(paste("\nFiltering proteins, this may take several minutes."))
filt_protein <- filtProt(irs_protein, # or spn_protein
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

# At this point there are no remaining missing values
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }


## ---- combine the final normalized data and sample meta data

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
	mutate(Subject = as.numeric(interaction(Mixture,Genotype))) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Subject, Intensity,Abundance) %>%
	group_by(Protein) %>%
	mutate(rel_Intensity = Intensity / sum(Intensity))


## ---- statistical analysis for overall Mutant-Controll contrast

# working with SL + IRS + SPN + (filt) + relative Intensity data
# ^ could be summarized as 4 key normalization steps + protein summarization (by sum,
# but could easily be changed to median) + peptide and protein level filtering

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores()-1)

# loop through all proteins, fit the model:
#fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)

# modeling subject as nested within mixture or seperately does not seem to have an effect on the results in this case
#fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Mixture:Subject) + (1|Subject)

## these seem to yield equivalent results...
fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Mixture:Subject)
#fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Subject)

# NOTE:
# adding the additional term SUBJECT changes our interpretation of the DF and variance
# fold change is the same

# example, fit the model to washc4 aka swip
swip <- gene_map$uniprot[grepl("Washc4", gene_map$symbol)]
fm <- lmerTest::lmer(fx, swip_tmt %>% subset(Protein == swip))

LT <- getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>%
	mutate(Contrast='Mutant-Control') %>%
	mutate(Protein = swip) %>%
        unique() %>% knitr::kable()


# loop to fit all proteins
proteins <- unique(swip_tmt$Protein)
results_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein == prot), control = lmer_control)
  LT <- getContrast(fm,"Mutant","Control")
  result = lmerTestContrast(fm,LT) %>%
  	     mutate(Contrast='Mutant-Control') %>%
	     mutate(Protein = prot) %>%
	     unique()
  return(result)
} #EOL


# collect results
df <- dplyr::bind_rows(results_list) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular) %>%
	arrange(Pvalue)


# check washc prot results
washc_prots <- gene_map$uniprot[grepl("Washc*",gene_map$symbol)]
df %>% filter(Protein %in% washc_prots) %>%
	mutate(FDR=formatC(FDR)) %>%
	mutate(Pvalue=formatC(Pvalue)) %>%
	dplyr::select(Protein, Symbol, Contrast, log2FC, SE, Tstatistic, Pvalue, FDR, DF) %>%
	knitr::kable()


# check nsig
sum(df$FDR<0.05)

head(df$Symbol)


# results for overall WT v MUT comparison
mut_wt_results <- df


## ---- FIT MODELS WITH SUBJECT TERM

swip_tmt <- swip_tmt %>% mutate(Subject = as.numeric(interaction(Mixture,Genotype)))

fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Subject)

fx <- log2(Intensity) ~ 1 + Condition + (1|Mixture) + (1|Mixture:Subject)


# loop
results_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein == prot), control = lmer_control)
  LT <- getContrast(fm,"Mutant","Control")
  result = lmerTestContrast(fm,LT) %>%
  	     mutate(Contrast='Mutant-Control') %>%
	     mutate(Protein = prot) %>%
	     unique()
  return(result)
} #EOL

# collect results
df <- dplyr::bind_rows(results_list) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular) %>%
	arrange(Pvalue)


# check nsig
sum(df$FDR<0.05)

# check top sig
head(df$Symbol)


## ---- what happens if SUBJECT is random?

swip_tmt <- swip_tmt %>% mutate(Subject = as.numeric(interaction(Mixture,Genotype)))

# randomize
swip_tmt$Subject <- sample(swip_tmt$Subject)

# loop
results_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein == prot), control = lmer_control)
  LT <- getContrast(fm,"Mutant","Control")
  result = lmerTestContrast(fm,LT) %>%
  	     mutate(Contrast='Mutant-Control') %>%
	     mutate(Protein = prot) %>%
	     unique()
  return(result)
} #EOL

# collect results
df <- dplyr::bind_rows(results_list) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular) %>%
	arrange(Pvalue)

# check nsig
sum(df$FDR<0.05)

# randomizing subject is equivalent to NOT including subject, makes sense


## ---- model SUBJECT nested within Mixture

swip_tmt <- swip_tmt %>% mutate(Subject = as.numeric(interaction(Mixture,Genotype)))

# nested model
fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Mixture:Subject)

# loop
results_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Protein == prot), control = lmer_control)
  LT <- getContrast(fm,"Mutant","Control")
  result = lmerTestContrast(fm,LT) %>%
  	     mutate(Contrast='Mutant-Control') %>%
	     mutate(Protein = prot) %>%
	     unique()
  return(result)
} #EOL

# collect results
df <- dplyr::bind_rows(results_list) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular) %>%
	arrange(Pvalue)

# check nsig
# nested Subject within Mixture
sum(df$FDR<0.05) # ~ not nested!



## ---- intra-BioFraction statistical analysis

doParallel::registerDoParallel(parallel::detectCores() - 1)

# all mut and wt conditions
biofractions <- c("F4","F5","F6","F7","F8","F9","F10")
mut <- paste("ConditionMutant",biofractions,sep=".")
wt <- paste("ConditionControl",biofractions,sep=".")

# loop
results_list <- foreach(prot = proteins) %dopar% {
  # fit the model to a given protein
  df <- swip_tmt %>%
	subset(Protein == prot) %>%
	mutate(rel_Intensity = Intensity/sum(Intensity))
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, df, control = lmer_control)
  # loop to evaluate n=7 contrasts for B=7 BioFractions
  n <- length(mut)
  contrasts <- list()
  res_list <- list()
  for (i in seq(n)) {
	LT <- getContrast(fm,mut[i],wt[i])
	res <- lmerTestContrast(fm,LT)
	res_list[[res$Contrast]] <- res
  }
  # collect results for all intra-BioFraction comparisons
  prot_results <- do.call(rbind, res_list) %>% mutate(Protein = prot)
  return(prot_results)
} #EOL


## ---- collect results for all intra-Biofraction comparisons

# collect results for all proteins
df <- dplyr::bind_rows(results_list) %>%
	group_by(Contrast) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular)

# collect as named list
results <- df %>% group_by(Contrast) %>% group_split()
names(results) <- sapply(results,function(x) unique(x$Contrast))

# shorten names
namen <- names(results)
shorter <- gsub("ConditionMutant\\.F[0-9]{1,2}-|ConditionControl\\.","",namen)
names(results) <- shorter

# sort and combine with overall Mutant-Control results
all_results <- results[biofractions]
class(all_results) <- "list"
all_results[["Mutant-Control"]] <- mut_wt_results

# sort
all_results <- lapply(all_results, function(x) x %>% arrange(Pvalue))

# summary of sig results
sapply(all_results,function(x) sum(x$FDR<FDR_alpha)) %>%
	t() %>% knitr::kable()

# all statistical results!
swip_results <-  do.call(rbind, all_results)

# sig_prots
temp_df <- swip_results %>% filter(Contrast=='Mutant-Control')
sig_prots <- unique(temp_df$Protein[temp_df$FDR<FDR_alpha])


## ----  Save key results

# save results as rda
myfile <- file.path(root,"data","swip_results.rda")
save(swip_results,file=myfile,version=2)
message("saved: ", myfile)

# write as excel
myfile <- file.path(root,"tables","SWIP-lmerTest-TMT-Results.xlsx")
write_excel(all_results, myfile)
message("saved: ", myfile)

# final normalized protein in tidy format as rda object
myfile <- file.path(datadir,"swip_tmt.rda")
save(swip_tmt,file=myfile,version=2)
message("saved: ", myfile)

# save gene_map
myfile <- file.path(datadir,"swip_gene_map.rda")
save(gene_map,file=myfile,version=2)
message("saved: ", myfile)

# save sig_prots
myfile <- file.path(datadir,"swip_sig_prots.rda")
save(sig_prots, file=myfile,version=2)
message("saved: ", myfile)
