#!/usr/bin/env Rscript 

# description: Working through MSstatsTMT protein-level models and 
# statistical comparisons

# input:
# * preprocessed protein-level data from PDtoMSstatsTMTFormat()
# comp - all intrafraction comparisons in the form "Mutant.F4-Control.F4"
# myfile <- file.path(root,"rdata","msstats_prot.rda")


## prepare environment --------------------------------------------------------

# project root dir and dir containing MSstatsTMT/R source code
root ="~/projects/SwipProteomics"
funcdir = "~/projects/SwipProteomics/src/MSstatsTMT"
#^ these are core MSstats internal functions used by MSstatsTMT_wrapper functions

# load renv
renv::load(root)

# load projects functions in root/R
# includes MSstatsTMT_wrappers.R utilized herein
devtools::load_all()

# load MSstatsTMT's guts
load_fun(funcdir)

# other imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
#	library(doParallel) # for parallel processing
})


# NOTE: attempts to parallelize fit and other functsion failed due to
# dependencies and OR fit. NOTE: even if all dependencies successfully passed
# to parallel processsing, then failure bc steps that perform statistics on fit
# needs to be done in same env that fit was created. This is due to the
# complexity of the underlying MSstatsTMT internal functions that do the 
# the statiscal operations on protein fits (calls to .updateModel).


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot


## build a contrast_matrix ----------------------------------------------------

# define all intrafraction comparisons:
comp <- paste(paste("Mutant",paste0("F",seq(4,10)),sep="."),
	      paste("Control",paste0("F",seq(4,10)), sep="."), sep="-")

# create a contrast matrix for given comparisons
conditions <- levels(msstats_prot$Condition)
contrast_matrix <- MSstatsTMT::getContrasts(comp, groups=conditions)

message("Contrast: ", rownames(contrast_matrix)[1])
knitr::kable(contrast_matrix[1,])


## begin protein-level modeling -------------------------------------------------

# for intra-fraction comparisons MSstatsTMT fits the model:
# FIXME: change run to batch or experiment... BATCH
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
message("Fitting protein-wise mixed-effects linear model of the form:\n\t",fx)

# Fit for Swip
prot <- "Q3UMB9"
fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, protein=prot)


## test contrasts -------------------------------------------------------------

# for each protein compare conditions declared in contrast_matrix
fit_list <- MSstatsTMT::testContrasts(fit_list, contrast_matrix, moderated = TRUE)

# get results list -- df for each comparison and compute padjust
results_list <- adjustPvalues(fit_list)

knitr::kable(bind_rows(results_list))


## loop to fit multiple proteins  -----------------------------------------------
# FIXME: implement pbar and error catching in fitLMER func

message("Fitting protein-wise mixed-effect linear models.")

#TODO: compile stats and compare to MSstats output
fit_list <- fitLMER(fx, msstats_prot, progress = TRUE)
# no applicable method for 'fixef' applied to an object of class "try-error"

fit_list <- testContrasts(fit_list, contrast_matrix, moderated = TRUE)
results_list <- adjustPvalues(fit_list)
knitr::kable(bind_rows(results_list))
