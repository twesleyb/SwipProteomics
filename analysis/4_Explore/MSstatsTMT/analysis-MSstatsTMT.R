#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT with MSstats
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# input data in root/rdata (too big for root/data as is)
input_psm = "5359_PSM_Report.xlsx"
input_samples = "5359_Samples.xlsx"

# Prepare the R environment ---------------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(MSstatsTMT)
	library(data.table)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })

# load sample data in root/rdata --------------------------------------

# pass colnames to read_excel
col_names <- c("S","Replicate","Sample","drop",
	       "Channel","ID","Condition","Mixture")
myfile <- file.path(root,"rdata",input_samples)
samples <- readxl::read_excel(myfile,col_names=col_names)

knitr::kable(head(samples))

# add colum for treatement 'Condition'
samples$Condition[grep("Control",samples$Condition)] <- "Control"
samples$Condition[grep("Mutant",samples$Condition)] <- "Mutant"
samples$Condition[grep("SPQC",samples$Condition)] <- "SPQC"

# drop un-needed col
samples$drop <- NULL


# read PSM data from excel ------------------------------------------------
myfile <- file.path(root,"rdata",input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)

# make columns look like MSstats by replacing special characters with '.'
chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
colnames(raw_pd) <- new_names

# add X if starts column starts with ".."
colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))

# build annotation file
length(unique(raw_pd$Spectrum.File)) # 36 unique spectrum files

length(unique(samples$Sample)) # there are 48 samples

## Map Spectrum.File to TMT Sample name ---------------------------------------

# annotation: data from which contains Run, Fraction, TechRepMixture, 
# Mixture, Channel, BioRepliicate, Condition (7) 
# Fraction: one technical replicate of one mixture may be fractionated into
# multiplee fractions to increasee analytical depth: 
# 1 tech rep of one mixture should cooresspond to multiple fractions Fraction =
# seq(1,12) = 12x fractions of first technical replicate cooresponding to 
# one biological subject, then they should have the same TechRepMixture and
# Mixture values.  Fraction =# 1,2,3,... 12

unique(raw_pd$Spectrum.File) # 36 unique spectrum files

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_", ID##### is experimental identifier.
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)

# collect all samples, grouped by Experiment.
all_samples <- samples$Sample
exp_samples <- split(all_samples, sapply(strsplit(all_samples,"_"),"[",1))

# There are 12 files per experiment.
x = sapply(exp_files,length)
knitr::kable(x)

# Each file cooresponds to the measurmeent of all 16 Fractions.
y = sapply(exp_samples,length)
knitr::kable(y)

# these twelve files are from experiment 1:
exp_files[[1]]

# and coorespond to these 16 samples:
exp_samples[[1]]

# Melt data
id_cols <- colnames(raw_pd)[!grepl("Abundance",colnames(raw_pd))]
tidy_pd <- reshape2::melt(raw_pd, id.vars = id_cols,
	       value.vars=colnames(raw_pd)[grepl("Abundance",colnames(raw_pd))])

# Complicated...
exp_to_samples <- lapply(names(exp_files), function(exp) {
	       x <- rep(paste(exp_samples[[exp]],collapse=";"),
			length(exp_files[[exp]]))
	       y <- setNames(x,nm=exp_files[[exp]])
	       return(y)
		    }) %>% unlist()

# 
names(exp_to_samples)
# exp_to_samples[1] == all 16 Samples delimited (;)

# map Spectrum.File to Samples
raw_pd$Spectrum.File <- exp_to_samples[raw_pd$Spectrum.File]

# split rows

# replace raw_pd$Spectrum.File with multiple samples
raw_pd$Spectrum.File

raw_pd$Spectrum.File



length(unique(knitr::kable(samples$Sample)))

# The format of Run in sample data does not match the pd data.
# Inspection shows that Run is interaction(experiment,Channel).
# According to MSstats docs, run should match input data 'Spectrum.File'
# Change raw_pd run annotatin to experiment.channel
exp <- sapply(strsplit(raw_pd$Spectrum.File,"_"),"[",1) # 3x experiments
fraction <- unlist(regmatches(raw_pd$Spectrum.File,
			     regexec("F[0-9]{1,2}",raw_pd$Spectrum.File)))

head(raw_pd$Spectrum.File)

#raw_pd$Spectrum.File <- interaction(exp,channel)


interaction(paste0("ID",samples$ID),samples$Channel)

annotation_pd <- data.table(Run="",
			    Fraction = "",
			    TechRepMixture = "",
			    Channel = "",
			    Condition = "",
			    Mixture = "",
			    BioReplicate = "")

annotation.pd$Run


# melt pd data into a tidy format
# 'Abundance #Channel' columns contain the numeric 'Intensity' data,
# all else are id.vars.
# What is this other 'Intensity' column?
colnames(raw_pd)[which(colnames(raw_pd)=="Intensity")] <- "PD.Intensity"
idy <- colnames(raw_pd)[!grepl("Abundance",colnames(raw_pd))]
tidy_pd <- reshape2::melt(raw_pd,id.vars=idy,
		  value.name="Intensity",variable.name="Channel")
tidy_pd$Channel <- gsub("Abundance: ","",tidy_pd$Channel)
# drop PSM's mapped to multiple proteins
tidy_pd <- tidy_pd %>% as.data.table() %>% filter(`# Proteins` == 1)
# collect required input for MSstats ----------------------------------
dt <- data.table(ProteinName = tidy_pd$"Master Protein Accessions",


# load the PSM data, coerce to MSstats format ---------------------------------

data_pd <- PDtoMSstats(psm_excel=file.path(root,"rdata",input_psm),
			samples_excel=file.path(root,"rdata",input_samples))

# The condition used for normalization must be "Norm"
data_pd$Condition[data_pd$Condition == "SPQC"] <- "Norm"


x = proteinSummarization(data_pd,
			 method="msstats",
			 global_norm=TRUE,
			 reference_norm=TRUE,
			 remove_norm_channel=TRUE)

data(raw.pd) # MSstats raw PD data
# keep psms with 1 matching protein
#psm_data[["filt"]] <- psm_data[["raw"]] %>% filter(`Number of Proteins` == 1)

# keep peptides with no missed cleavage
#psm_data[["no missed cleavage"]] <- psm_data[["filt"]] %>% filter(`Number of Missed Cleavages` == 0)

# summary:
#knitr::kable(t(formatC(sapply(psm_data,nrow),big.mark=",")))

#data(input.pd) # MSstatsTMT dataset

# tidy_peptide from root/rdata, its big ~ 99 mb
#load(file.path(root,"rdata","tidy_peptide.rda"))

#str(tidy_peptide)
#str(input.pd)

## build input to MSstats -----------------------------------------------------
## REALLY NEED MORE INFO ABOUT THE COLUMNS
# the dta should contain the following columns:
# * ProteinName - uniprot Accession
# * PSM - peptide spectrum match = an ionized peptide, what is measured
# * TechRepMixture - QC?
# * Mixture
# * Run
# * Channel
# * Condition
# * BioReplicate
# * Intensity

# Number of Missed Cleavages --> remove peptides with missed cleavage?
# Number of Proteins --> remove multiples?
# total Ions
# ions matched
#raw_psm$"Number of Proteins"[94]
#raw_psm$"Master Protein Accessions"[94]
#all(sapply(psm_data[["no missed cleavage"]][["Master Protein Accessions"]],length)==1)

df <- psm_data[["no missed cleavage"]]
# exctact sample info from "Spectrum File"
ids <- gsub("ID","", sapply(strsplit(df$"Spectrum File","_"),"[",1))
j
samples$"Proteomics ID" <- sapply(strsplit(samples$Sample,"\\."),"[",5)

unique(df[["Spectrum File"]])

input_df <- data.table("ProteinName" = as.factor(df[["Master Protein Accessions"]]), # Uniprot Accession
		 "PSM" = as.factor(df[["Annotated Sequence"]]),
		 "TechRepMixture" = as.factor(tidy_peptide[["Fraction"]]),
		 "Mixture" = as.factor(tidy_peptide[["Experiment"]]), # Exp1, Exp2, Exp3
		 "Run" = as.factor(tidy_peptide[["Sample"]]), # Sample Name
		 "Channel" = as.factor(tidy_peptide[["Channel"]]), # 126, 127C, ect
		 "Condition" = as.factor(tidy_peptide[["Treatment"]]), # SPQC, WT, MUT
		 "BioReplicate" = as.factor(paste(tidy_peptide[["Experiment"]],tidy_peptide[["Treatment"]],sep="_")),
		 "Intensity" = as.numeric(tidy_peptide[["Intensity"]]))
