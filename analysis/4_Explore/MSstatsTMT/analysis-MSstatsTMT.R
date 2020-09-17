#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT with MSstats
# author: Tyler W Bradshaw <twesleyb10@gmail.com>


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
#root <- Sys.getEnv("ROOT") # set by activate script in 
#FIXME: where is 'Error in readRDS(file) : unknown input format' coming from?
root <- getrd()
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(MSstats)
	library(MSstatsTMT)
	library(data.table)
})

# load functions in root/R 
suppressWarnings({ devtools::load_all() })

# load data in root/data
myfile <- list.files(file.path(root,"data"),pattern="PSM",full.names=TRUE)
psm_data <- list("raw"=data.table::fread(myfile))

# keep psms with 1 matching protein
psm_data[["filt"]] <- psm_data[["raw"]] %>% filter(`Number of Proteins` == 1)

# keep peptides with no missed cleavage
psm_data[["no missed cleavage"]] <- psm_data[["filt"]] %>% filter(`Number of Missed Cleavages` == 0)

# summary:
knitr::kable(t(formatC(sapply(psm_data,nrow),big.mark=",")))

#data(input.pd) # MSstatsTMT dataset

# tidy_peptide from root/rdata, its big ~ 99 mb
#load(file.path(root,"rdata","tidy_peptide.rda")) 

#str(tidy_peptide)
#str(input.pd)

## build input to MSstats -----------------------------------------------------
## REALLY NEED MORE INFO ABOUT THE COLUMNS
# the data should contain the following columns:  
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
