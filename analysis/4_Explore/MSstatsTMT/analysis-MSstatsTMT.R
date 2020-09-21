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

# load the sample data ---------------------------------------------------

# load sample data in root/rdata
myfile <- file.path(root,"rdata",input_samples)
col_names <- c("S","Replicate","Run","Sample",
	       "Channel","ID","Condition","Mixture")
samples <- readxl::read_excel(myfile,col_names=col_names)

# clean-up sample info
samples$S <- NULL
samples$Sample <- NULL
samples$Replicate <- gsub("F","",samples$Replicate)
samples$Condition[grep("Control",samples$Condition)] <- "Control"
samples$Condition[grep("Mutant",samples$Condition)] <- "Mutant"
samples$Condition[grep("SPQC",samples$Condition)] <- "SPQC"
samples$BioReplicate <- interaction(samples$Condition,
				    samples$Replicate)

# build annotation df:
# Run, Fraction, TechRepMixturee, Channel, Conditon, Bioreplicate,Mixture
anno_df <- data.table(Run=as.factor(samples$ID),
		      Fraction=rep(1,nrow(samples)),
		      TechRepMixture=as.factor(samples$Channel),
		      Channel=as.factor(samples$Channel),
		      Condition=as.factor(samples$Condition),
		      BioReplicate=as.factor(samples$BioReplicate),
		      Mixture=as.factor(samples$Mixture))

# load the PSM data ---------------------------------------------------

# load PSM data in root/rdata
myfile <- file.path(root,"rdata",input_psm)
psm_df <- readxl::read_excel(myfile)

x = psm_data[[1]]

# coerce data to MSstats format -----------------------------------------------

df <- PDtoMSstatsTMTFormat(psm_df, anno_df,
		     rmPSM_withMissing_withinRun=TRUE)

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


x = proteinSummarization(input_df,
			 method="msstats",
			 global_norm=TRUE,
			 reference_norm=TRUE,
			 remove_norm_channel=TRUE)
