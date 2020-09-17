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
	library(MSstats)
	library(MSstatsTMT)
	library(data.table)
	library(dplyr)
})

# load functions in root/R 
suppressWarnings({ devtools::load_all() })

# load data in root/data
data(input.pd) # MSstatsTMT dataset

# tidy_peptide from root/rdata, its big ~ 99 mb
load(file.path(root,"rdata","tidy_peptide.rda")) 


str(tidy_peptide)
str(input.pd)

quit()

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

input_df <- data.table("Protein" = as.factor(tidy_peptide[["Accession"]]), # Uniprot Accession
		 "PSM" = as.factor(paste(tidy_peptide[["Sequence"]],tidy_peptide[["Modifications"]],sep="_")),
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
