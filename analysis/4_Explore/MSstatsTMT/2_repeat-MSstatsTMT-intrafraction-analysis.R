#!/usr/bin/env Rscript 

# title:
# author: twab
# description: Working through MSstatsTMT protein-level statistics

# input:


# options:
n = 100
save_rda = FALSE


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
load(file.path(root,"rdata","msstats_prot.rda"))
load(file.path(root,"rdata","msstats_contrasts.rda"))


## begin protein-level modeling -------------------------------------------------
# FIXME: change run to batch or experiment... BATCH

# for intra-fraction comparisons MSstatsTMT fits the model:
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
message("Fitting protein-wise mixed-effects linear model of the form:\n\t",fx)

# | Run cooresponds to a TMT multiplex experiment aka a Batch or Experiment.

# subset the data
proteins <- unique(as.character(msstats_prot$Protein))
if (n > length(proteins) | n == 'all' ) {
	# no subset
	message("No subset.")
	prots = proteins
	message("Analyzing all (N=",length(proteins),") proteins.")
} else {
	prots = sample(proteins,n)
	message("Analyzing subset of proteins (n=",n,").")
}

# start timer here:
t0 = Sys.time()

# fit the protein-wise models
fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, protein=prots)

## test contrasts -------------------------------------------------------------

# for each protein compare conditions declared in contrast_matrix
fit_list <- MSstatsTMT::testContrasts(fit_list, msstats_contrasts, moderated = TRUE)

# get results list -- df for each comparison and compute padjust
results_list <- MSstatsTMT::getResults(fit_list)
t1 = Sys.time()

message(difftime(t1-t0,"seconds"))


quit()



## loop to fit multiple proteins  -----------------------------------------------
#TODO: compile stats and compare to MSstats output

message("Fitting protein-wise mixed-effect linear models.")

# fit models
# FIXME: how to improve speed?
#Model failed to converge with max|grad| = 0.00223701 (tol = 0.002, component 1) 
fit_list <- fitLMER(fx, msstats_prot, progress = TRUE)

fit_list <- testContrasts(fit_list, contrast_matrix, 
			  moderated = TRUE, progress = TRUE)

# test contrasts
#69%Error in Lc %*% as.matrix(vcov_out) : non-conformable arguments 
# no error when looping like this???:
# no errors this time...
for (i in c(1:length(fit_list))) {
	print(i)
fits[[i]] <- testContrasts(fit_list[i], contrast_matrix, 
			  moderated = TRUE, progress = FALSE)
}

fit_list <- unlist(fits, recursive=FALSE)

# compile results
results_list <- adjustPvalues(fit_list)

# summary of significant proteins:
sapply(results_list,function(x) sum(x$"P-adjust" < 0.05,na.rm=TRUE))

## ----------------------------------------------------------------------------
# save
myfile <- file.path(root,"rdata","msstats_intrafraction_results.rda")
save(results_list,file=myfile,version=2)
