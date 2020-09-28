#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT with MSstats
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# input data in root/rdata (too big for root/data as is)
#input_dir = "5359_PSM_Data.zip"
input_psm = "rdata/5359_PSM_Report.xlsx"
input_samples = "rdata/5359_Sample_Report.xlsx"


# Prepare the R environment ---------------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
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

# FIXME: todo


# load sample data in root/rdata --------------------------------------

# this excel spreadsheet was from Greg, exported from PD(?)
myfile <- file.path(root,input_samples)

# pass colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
samples <- readxl::read_excel(myfile,col_names=col_names)

## read PSM data from excel ------------------------------------------------

# read raw data, NOTE: this takes several minutes
myfile <- file.path(root,input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)


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
samples$ConditionFraction <- NULL

# Set the formatting of the channel for normalization for MSstats
samples$Condition[samples$Condition=="SPQC"] <- "Norm" 

# clean-up Mixtrure column 
samples$Mixture <- gsub("F","M",samples$Mixture)

# Remove un-needed col
samples$drop <- NULL

# FIXME: how should BioReplicate be defined?
samples$BioReplicate <- paste0(samples$Condition,
			       as.numeric(as.factor(samples$Experiment)))


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
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$Fraction <- named_fractions[annotation_dt$MS.Channel]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
annotation_dt$Subject <- samples$BioReplicate[idx]
annotation_dt$Channel <- samples$Channel[idx]

# how to handle repeated measures design? Subject.Fraction
annotation_dt$BioReplicate <- as.character(interaction(samples$BioReplicate[idx],
						       samples$Fraction[idx]))

# whoops, fix norm#.F# -- should be just norm - the spqc sample
annotation_dt$BioReplicate[grepl("norm",annotation_dt$BioReplicate)] <- "Norm"

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL

# check the annotation file
knitr::kable(annotation_dt[c(1:16),])

# save to file
myfile <- file.path(root,"rdata","annotation_dt.rda")
save(annotation_dt,file=myfile,version=2)
message(paste("saved",myfile))

## coerce data to MSstats format ---------------------------------

# NOTE: this takes a considerable amount of time
data_pd <- PDtoMSstatsTMTFormat(raw_pd,annotation_dt)

# ** Shared PSMs (assigned in multiple proteins) are removed.
# ** 705 features have 1 or 2 intensities across runs and are removed.
# ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.
# ** For peptides overlapped between fractions of M1_1, use the fraction with maximal average abundance.
# ** For peptides overlapped between fractions of M2_1, use the fraction with maximal average abundance.
# ** For peptides overlapped between fractions of M3_1, use the fraction with maximal average abundance.
# ** Fractions belonging to same mixture have been combined.

#knitr::kable(head(data_pd))

# save to file
myfile <- file.path(root,"rdata","data_pd.rda")
save(data_pd,file=myfile,version=2)
message(paste("saved",myfile))


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
# Normalization between MS runs for Protein : Q9DA08 ( 8525  of  8556 )
data_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)

# ** Protein-level summarization done by MSstats.
# There were 14 warnings
# 1: In survreg.fit(X, Y, weights, offset, init = init, controlvals = control...
# Ran out of iterations and did not converge
# Please check it. At this moment, normalization is not performed.

# save to file
myfile <- file.path(root,"rdata","data_prot.rda")
save(data_prot,file=myfile,version=2)
message(paste("saved",myfile))
	

# Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
all_results <- groupComparisonTMT(data_prot)	

# get boundary fit error from  lme... something is not correct with our design
# or MSstats's interpretation of the design.
