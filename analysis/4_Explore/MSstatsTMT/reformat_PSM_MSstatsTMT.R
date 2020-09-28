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

# this excel spreadsheet was from Greg, exported from PD?
# pass colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
myfile <- file.path(root,"rdata",input_samples)
samples <- readxl::read_excel(myfile,col_names=col_names)

# Extract relevant annotations for MSstatsTMT
# Format colums for MSstatsTMT
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

# munge extra

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
annotation_dt$BioReplicate <- as.character(interaction(samples$BioReplicate[idx],samples$Fraction[idx]))

# whoops, fix norm#.F# -- should be just norm - the spqc sample
annotation_dt$BioReplicate[grepl("norm",annotation_dt$BioReplicate)] <- "norm"

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL

# knitr::kable(annotation_dt[c(1:16),])
# |Run                                                              | Fraction| TechRepMixture|Mixture |Condition |Subject  |Channel |BioReplicate |
# |:----------------------------------------------------------------|--------:|--------------:|:-------|:---------|:--------|:-------|:------------|
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        1|              1|M1      |Control   |Control1 |126     |Control1.F4  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        3|              1|M1      |Mutant    |Mutant1  |127N    |Mutant1.F4   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        2|              1|M1      |Control   |Control1 |127C    |Control1.F5  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        5|              1|M1      |Mutant    |Mutant1  |128N    |Mutant1.F5   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        4|              1|M1      |Control   |Control1 |128C    |Control1.F6  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        7|              1|M1      |Mutant    |Mutant1  |129N    |Mutant1.F6   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        6|              1|M1      |Control   |Control1 |129C    |Control1.F7  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        9|              1|M1      |Mutant    |Mutant1  |130N    |Mutant1.F7   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |        8|              1|M1      |Control   |Control1 |130C    |Control1.F8  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       11|              1|M1      |Mutant    |Mutant1  |131N    |Mutant1.F8   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       10|              1|M1      |Control   |Control1 |131C    |Control1.F9  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       13|              1|M1      |Mutant    |Mutant1  |132N    |Mutant1.F9   |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       12|              1|M1      |Control   |Control1 |132C    |Control1.F10 |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       15|              1|M1      |Mutant    |Mutant1  |133N    |Mutant1.F10  |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       14|              1|M1      |norm      |norm1    |133C    |norm         |
# |ID53395_01_FAIMS_OTIT_FSN20357_5539_021120_F11.raw               |       16|              1|M1      |norm      |norm1    |134N    |norm         |

# save to file
myfile <- file.path(root,"rdata","annotation_dt")
save(annotation_dt,file=myfile,version=2)
message(paste("saved",myfile))

## coerce data to MSstats format ---------------------------------

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

# FIXME: Joining, by = (Run, Channel )
data_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE,	
                                  remove_empty_channel = TRUE)	


# save to file
myfile <- file.path(root,"rdata","data_prot")
save(data_prot,file=myfile,version=2)
message(paste("saved",myfile))
	

# Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
all_results <- groupComparisonTMT(data_prot)	

# The linear models utilized by MSstats:
# [1]



