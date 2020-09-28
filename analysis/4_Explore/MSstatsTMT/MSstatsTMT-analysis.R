#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------

# input data in root/data:
input_dir = "data/PSM.zip"

# PSM.zip contains:
input_psm = "PSM/5359_PSM_Report.xlsx"
input_samples = "PSM/5359_Sample_Report.xlsx"

## Prepare the working environment --------------------------------------------

root <- "~/projects/SwipProteomics"
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata") # temporary/large files not tracked by git

# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(MSstatsTMT)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## unzip the directory containing the raw data ---------------------------------

unzip(file.path(root,input_dir)) # unzip root/PSM.zip into current dir


# load sample data in root/rdata --------------------------------------
# this excel spreadsheet was from Greg, exported from PD(?)

# pass colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
samples <- readxl::read_excel(input_samples,col_names=col_names)

unlink(input_samples)


## read PSM data from excel ------------------------------------------------

# read raw data, NOTE: this takes several minutes
message("\nLoading raw PSM data.") 
raw_pd <- readxl::read_excel(input_psm,progress=FALSE)

unlink(input_psm) 

unlink(tools::file_path_sans_ext(basename(input_dir))) # rmdir ./PSM

## re-format sample annotations for MSstats -----------------------------------
# Extract relevant annotations for MSstatsTMT
# Format colums for MSstatsTMT

f1 <- function(x) { # x = samples$ConditionFraction
	# extract sample 'Condition'
	paste0("F",as.numeric(sapply(strsplit(x,"Control|Mutant|SPQC"),"[",2)))
}

f2 <- function(x) { # x = samples$ConditionFraction 
	# extract sample 'Fraction'
	sapply(strsplit(x,"[0-9]{1,2}"),"[",1)
}

# Munge ConditionFraction Column to Fraction and Condition cols
samples$Fraction <- f1(samples$ConditionFraction)
samples$Condition <- f2(samples$ConditionFraction)
samples$Condition <- as.character(interaction(samples$Condition,samples$Fraction)) # added
samples$ConditionFraction <- NULL

# Set the formatting of the channel for normalization for MSstats
samples$Condition[grepl("SPQC",samples$Condition)] <- "Norm" 

# clean-up Mixture column 
samples$Mixture <- gsub("F","M",samples$Mixture)

# Remove un-needed col
samples$drop <- NULL

# FIXME: how should BioReplicate be defined?
samples$BioReplicate <- paste(samples$Condition,paste0("R",as.numeric(as.factor(samples$Experiment))),sep=".")


## re-format PSM data for MSstatsTMT ---------------------------------------------

# make columns look like MSstats by replacing special characters with '.'
chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
colnames(raw_pd) <- new_names

# add X if starts column starts with ".."
colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))


## map Spectrum.Files to MS.Channels ------------------------------------------

# munge extra

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_", ID##### is experimental identifier.
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment ID" = rep(names(exp_files),
					     times=sapply(exp_files,length)),
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

# annotation: data.frame which contains 
# Run - match Spectrum.File
# Fraction - TMT mixture may be fractionated into multiple fractions
# TechRepMixture - all 1, no repeats of a mixture
# Mixture - concatenation of TMT labeled samples - an MS experiment
# Channel - the TMT channels
# BioReplicate - individual subjects (mice)
# Condition - WT MUT norm (SPQC)

# create annotation_dt from Spectrum.Files and MS.Runs
message("\nBuilding annotation data.frame for MSstats.") 
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$Fraction <- named_fractions[annotation_dt$MS.Channel]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
#annotation_dt$Subject <- samples$BioReplicate[idx]
annotation_dt$Channel <- samples$Channel[idx]

# Add BioFraction column 
# FIXME: how to pass additional covariates to MSstats?
annotation_dt$BioFraction <- samples$Fraction[idx]

# how to handle repeated measures design? 
annotation_dt$BioReplicate <- samples$BioReplicate[idx]
#as.character(interaction(samples$BioReplicate[idx],
#						       samples$Fraction[idx]))

# whoops, fix norm#.F# -- should be just Norm - the SPQC sample which is
# analyzed in Technical duplicate for every experiment?
annotation_dt$BioReplicate[grepl("Norm",annotation_dt$BioReplicate)] <- "Norm"

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL

# check the annotation file
# this basic design is repeated 3x (M = 3):
knitr::kable(annotation_dt[c(1:16),])

# save to file
myfile <- file.path(datadir,"annotation_dt.rda")
save(annotation_dt,file=myfile,version=2)
message(paste("\nSaved",myfile))


## coerce data to MSstats format ---------------------------------

# NOTE: this takes a considerable amount of time
message("\nConverting PSM data to MSstatsTMT format.")
data_pd <- PDtoMSstatsTMTFormat(raw_pd, annotation_dt)

# save to file
myfile <- file.path(rdatdir,"data_pd.rda")
save(data_pd,file=myfile,version=2)
message(paste("\nSaved",myfile))


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
message("\nPerforming normalization and protein sumamrization.")
data_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)

# save to file
myfile <- file.path(rdatdir,"data_prot.rda")
save(data_prot,file=myfile,version=2)
message(paste("\nSaved",myfile))
	
quit()

# Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
all_results <- groupComparisonTMT(data_prot)	

myfile <- file.path(rdatdir,"all_results.rda")
save(all_results,file=myfile,version=2)
message(paste("\nSaved",myfile))

fwrite(all_results,"MSstatsTMT_Results.csv")
