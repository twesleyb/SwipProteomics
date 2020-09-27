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

# this is data from Greg, exported from PD?
# pass colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
myfile <- file.path(root,"rdata",input_samples)
samples <- readxl::read_excel(myfile,col_names=col_names)
samples$Fraction <- paste0("F",as.numeric(sapply(strsplit(samples$ConditionFraction,"Control|Mutant|SPQC"),"[",2)))
samples$Condition <- sapply(strsplit(samples$ConditionFraction,"[0-9]{1,2}"),"[",1)
samples$ConditionFraction <- NULL
samples$Condition[samples$Condition=="SPQC"] <- "norm"
samples$Mixture <- gsub("F","M",samples$Mixture)
samples$drop <- NULL
samples$BioReplicate <- paste0(samples$Condition,as.numeric(as.factor(samples$Experiment)))


# read PSM data from excel ------------------------------------------------

# read raw data, this takes several minutes
myfile <- file.path(root,"rdata",input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)

# make columns look like MSstats by replacing special characters with '.'
chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
colnames(raw_pd) <- new_names

# add X if starts column starts with ".."
colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))


## map Spectrum.Files to MS.Channels ------------------------------------------

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_", ID##### is experimental identifier.
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment ID" = rep(names(exp_files),times=sapply(exp_files,length)),
	   "Run" = unlist(exp_files))

# collect all MS.Channels, grouped by Experiment
all_channels <- samples$MS.Channel
exp_channels <- split(all_channels, sapply(strsplit(all_channels,"_"),"[",1))
exp_dt <- data.table("Experiment ID" = rep(names(exp_channels),
					times=sapply(exp_channels,length)),
	   "MS.Channel" = unlist(exp_channels))

# use exp_channels to create numeric ID for MS Fraction
exp_fraction_list <- lapply(exp_channels, function(x) {
	       setNames(as.numeric(as.factor(x)),x)
	   })
x <- unlist(exp_fraction_list)
named_fractions <- setNames(x,sapply(strsplit(names(x),"\\."),"[",2))


# build annotation file -------------------------------------------------------

# annotation: data from which contains 
# Run - match Spectrum.File
# Fraction - TMT mixture may be fractionated into multiple fractions
# TechRepMixture - all 1, no repeats of a mixture
# Mixture - concatenation of TMT labeled samples - an MS experiment
# Channel - the TMT channels
# BioReplicate - 
# Condition (7)  - WT MUT norm (SPQC)
## Add BioFraction - the subcellular fraction
# BioFraction - 
# Combine (spectrum.)files_dt with exp_dt containing MS.Runs 

annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$Fraction <- named_fractions[annotation_dt$MS.Channel]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
annotation_dt$Channel <- samples$Channel[idx]
annotation_dt$BioReplicate <- samples$BioReplicate[idx]
annotation_dt$BioFraction <- samples$Fraction[idx]

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL
annotation_dt$"BioFraction" <- NULL


## coerce data to MSstats format ---------------------------------

data_pd <- PDtoMSstatsTMTFormat(raw_pd,annotation_dt)
