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

#length(unique(samples$Sample)) # there are 48 samples

## Map Spectrum.File to TMT Sample name ---------------------------------------

# annotation: data from which contains Run, Fraction, TechRepMixture, 
# Mixture, Channel, BioRepliicate, Condition (7) 
# Fraction: one technical replicate of one mixture may be fractionated into
# multiplee fractions to increasee analytical depth: 
# 1 tech rep of one mixture should cooresspond to multiple fractions Fraction =
# seq(1,12) = 12x fractions of first technical replicate cooresponding to 
# one biological subject, then they should have the same TechRepMixture and
# Mixture values.  Fraction =# 1,2,3,... 12

# spectrum.file channel 126, KO
# spectrum.file channle 127, WT

#unique(raw_pd$Spectrum.File) # 36 unique spectrum files

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_", ID##### is experimental identifier.
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment" = rep(names(exp_files),times=sapply(exp_files,length)),
	   "Spectrum.File" = unlist(exp_files))

# collect all samples, grouped by Experiment.
all_samples <- samples$Sample
exp_samples <- split(all_samples, sapply(strsplit(all_samples,"_"),"[",1))
exp_dt <- data.table("Experiment" = rep(names(exp_samples),times=sapply(exp_samples,length)),
	   "MS.Run" = unlist(exp_samples))
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment")

# there should be M mixtures x F fractions x C channels = 3 x 12 x 16 = 576 rows

# annotate
idx <- match(annotation_dt$MS.Run, samples$Sample) # FIXME: Rename to MS.Run or MS.Channel

annotation_dt


# load the PSM data, coerce to MSstats format ---------------------------------

data_pd <- PDtoMSstats(psm_excel=file.path(root,"rdata",input_psm),
			samples_excel=file.path(root,"rdata",input_samples))

