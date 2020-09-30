#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
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
	library(data.table)
})

# load functions in root/R
suppressWarnings({ devtools::load_all() })

# load sample data in root/rdata --------------------------------------

# pass colnames to read_excel
col_names <- c("S","Replicate","Run","Sample",
	       "Channel","ID","Condition","Mixture")
myfile <- file.path(root,"rdata",input_samples)
samples <- readxl::read_excel(myfile,col_names=col_names)

# clean-up sample info treatement 'Condition'
samples$Condition[grep("Control",samples$Condition)] <- "Control"
samples$Condition[grep("Mutant",samples$Condition)] <- "Mutant"
samples$Condition[grep("SPQC",samples$Condition)] <- "SPQC"

# read PSM data from excel ------------------------------------------------
myfile <- file.path(root,"rdata",input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)


