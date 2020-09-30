#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT with MSstats
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# input data in root/rdata (too big for root/data as is)
input_psm = "5359_PSM_Report.xlsx"
input_samples = "5359_Samples.xlsx"

## Misc Functions -------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# get the repository's root directory
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

psm_file = file.path(root,"rdata",input_psm)
samples_file = file.path(root,"rdata",input_samples)
data_pd <- PDtoMSstats(psm_file,samples_file)

PDtoMSstats <- function(psm_file,samples_file) {
	# a function that coerces PSM data to MSstats format
	# input is path to excel psm report and excel samples report
	# FIXME: need to drop ambiguous, drop rows with missing abundance,
	# summarize multiple psms
	# load sample data in root/rdata --------------------------------------
	col_names <- c("S","Replicate","Run","Sample",
		       "Channel","ID","Condition","Mixture")
	samples <- readxl::read_excel(samples_file,col_names=col_names)
	# clean-up sample info treatement 'Condition'
	samples$Condition[grep("Control",samples$Condition)] <- "Control"
	samples$Condition[grep("Mutant",samples$Condition)] <- "Mutant"
	samples$Condition[grep("SPQC",samples$Condition)] <- "SPQC"
	# read PSM data from excel and melt into tidy format -------------------
	raw_pd <- readxl::read_excel(psm_file,progress=FALSE)
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
	        PeptideSequence = tidy_pd$"Sequence",
	        Charge = tidy_pd$"Charge",
	        PSM = paste(tidy_pd$Sequence,tidy_pd$Charge,sep="_"),
	        Mixture = sapply(strsplit(tidy_pd$"Spectrum File","_"),"[",1),
	        TechRepMixture=1,
	        Run = "",
	        Channel=tidy_pd$"Channel",
	        BioReplicate=gsub("\\.",".F",gsub("F","Exp",tidy_pd$"File ID")),
	        Intensity=tidy_pd$"Intensity",
	        Condition="")
	# add Run annotation as Mixture(experiment/ID).Channel
	dt$Run <- paste(dt$Mixture,dt$Channel,sep=".")
	# use 'samples' to annotate PD data with treatment 'Condition'
	samples$Mixture.Channel <- paste(sapply(strsplit(samples$Run,"_"),"[",1),
					 samples$Channel,sep=".")
	idx <- match(dt$Run,samples$Mixture.Channel)
	dt$Condition <- samples$Condition[idx]
	# summarize PSMs with multiple measurements
	#foo = dt %>% group_by(Run,PSM) %>% summarize(Intensity = Sum(Intensity))
	# return the data reformated for analysis with MSstats ----------------
	return(dt)
}


# Prepare the R environment ---------------------------------------------------

# load renv
root <- getrd()
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(MSstatsTMT)
	library(data.table)
})

# load functions in root/R
suppressWarnings({ devtools::load_all() })

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
