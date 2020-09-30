#!/usr/bin/env Rscript

#' ---	
#' title: MSstatsTMT 
#' A package for analysis of significance in shotgun mass spectrometry-based
#' proteomic experiments with tandem mass tag (TMT) labeling"	
#' author: Ting Huang, Meena Choi, Sicheng Hao, Olga Vitek
#' ---	

#' This vignette summarizes the introduction and various options of all
#' functionalities in MSstatsTMT. 	

#' MSstatsTMT includes the following three steps for statistical testing: 	

#' 1. Converters for different peptide quantification tools to get the input
#' with required format: `PDtoMSstatsTMTFormat`, `MaxQtoMSstatsTMTFormat`,
#' `SpectroMinetoMSstatsTMTFormat` and `OpenMStoMSstatsTMTFormat`.	
#' 2. Protein summarization based on peptide quantification data:
#' `proteinSummarization`	
#' 3. Group comparison on protein quantification data:  `groupComparisonTMT`	

## functions ------------------------------------------------------------------

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


knit_str <- function(robject) {
	# a helper function to remove repeating whitespace
	removews <- function(string){
		return(gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE))
	}
	# caputure output of str
	raw <- capture.output(str(robject))
	# init a data.frame that will contain info summarizing the robject
	df <- data.frame("Summary" = gsub("\t", " ",raw[1])) # remove tabs
	# parse variable names
	var_namen <- trimws({
		gsub("\\$","",sapply(strsplit(raw[2:length(raw)],":"),"[",1)) })
	# parse descriptions of the variables
	vars <- trimws({
		gsub("\\$","",sapply(strsplit(raw[2:length(raw)],":"),"[",2)) })
	# remove trailing ... and truncate to 60 characters
	clean_vars <- substr(trimws(gsub("\\...|,.."," ",vars)),1,60)
	res <- cbind(df,data.frame(setNames(as.list(clean_vars),nm=var_namen)))
	# pretty print the result
	knitr::kable(t(res),col.names="")
}

# prepare the working environment ----------------------------------------------
# load renv
root <- getrd()
renv::load(root)

# imports
library(MSstatsTMT)	

# load the data
data(raw.pd)
data(annotation.pd)


# PDtoMSstatsTMTFormat --------------------------------------------------------
# Preprocess PSM data from Proteome Discoverer and convert into the required
# input format for MSstatsTMT.	
	
# coerce input to MSstats format
input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)	


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
quant.msstats <- proteinSummarization(input.pd,	
                                      method="msstats",	
                                      global_norm=TRUE,	
                                      reference_norm=TRUE,	
                                      remove_norm_channel = TRUE,	
                                      remove_empty_channel = TRUE)	
	

# Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
test.pairwise <- groupComparisonTMT(quant.msstats)	
	
# Only compare condition 0.125 and 1	
comparison<-matrix(c(-1,0,0,1),nrow=1)	
# Set the names of each row	
row.names(comparison)<-"1-0.125"	
# Set the column names	
colnames(comparison)<- c("0.125", "0.5", "0.667", "1")	
	
test.contrast <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison)	
