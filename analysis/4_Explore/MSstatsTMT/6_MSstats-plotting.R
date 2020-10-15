#!/usr/bin/env Rscript

# title:
# author: twab
# description: 

# load renv
root = "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(MSstats)
  library(MSstatsTMT)
})

# local imports
devtools::load_all(quiet=TRUE)

# load the data ---------------------------------------------------------------

data(swip)

myfile<-file.path(root,"rdata","msstats_input.rda")
data(
load(myfile) 

myfile<-file.path(root,"data","msstats_prot.rda")
load(myfile) 

## dataProcessPlotsTMT() ------------------------------------------------------	
# (1) profile plot (specify "ProfilePlot" in option type), to identify the
# potential sources of variation for each protein;	
# (2) quality control plot (specify "QCPlot" in option type), to evaluate the
# systematic bias between MS runs.	

## Profile plot with all the channels	
## Profile plot without norm channnels and empty channels	
dataProcessPlotsTMT(data.peptide = msstats_input,	
                     data.summarization = msstats_prot,	
                     type = 'ProfilePlot', # or QCPlot
		     which.Protein = swip, #"all", # or allonly for QCPlot
		     originalPlot = TRUE,
		     summaryPlot = TRUE, 
                     width = 42, # adjust the figure width since there are 15 TMT runs.	
                     height = 7,
		     address=file.path(root,"figs","MSstatsTMT"))	
	
## Quality control plot 	
dataProcessPlotsTMT(data.peptide = msstats_input,	
                     data.summarization = msstats_prot,	
                     type = 'QCPlot', # or QCPlot
		     which.Protein = swip, #"all", # or allonly for QCPlot
		     originalPlot = TRUE,
		     summaryPlot = TRUE, 
                     width = 42, # adjust the figure width since there are 15 TMT runs.	
                     height = 7,
		     address=file.path(root,"figs","MSstatsTMT"))	

## NOTE:
# * originalPlot = TRUE(default): draw original profile plots, 
#   without normalization
# * summaryPlot = TRUE(default): draws profile plots with protein summarization
#   for each channel and MS run.
