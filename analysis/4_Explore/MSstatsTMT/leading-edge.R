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
	library(doParallel) # for parallel processing
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot


## build a contrast_matrix ----------------------------------------------------

# define all intrafraction comparisons:
comp <- paste(paste("Mutant",paste0("F",seq(4,10)),sep="."),
	      paste("Control",paste0("F",seq(4,10)), sep="."), sep="-")

message("Intrafraction Comparisons:")
knitr::kable(comp)

# create a contrast matrix for given comparisons
contrast_matrix <- getContrasts(comp, groups=levels(msstats_prot$Condition))

dim(contrast_matrix)
rownames(contrast_matrix)[1]
knitr::kable(contrast_matrix[1,])


## begin protein-level modeling -------------------------------------------------

# for intra-fraction comparisons MSstatsTMT fits the model:
# FIXME: change run to batch or experiment... BATCH
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
message("Fitting protein-wise mixed-effects linear model of the form:\n\t",fx)

# fit model for 9x + Swip
# FIXME: timit! parallize!
#prots <- c("Q3UMB9", sample(unique(as.character(msstats_prot$Protein)),9))
proteins <- unique(as.character(msstats_prot$Protein))
proteins <- sample(proteins,100)

#fit_list <- fitLMER(fx,msstats_prot,protein=prots)

# Register some nodes to do work.
nThreads <- parallel::detectCores() - 1
workers <- parallel::makeCluster(nThreads, type = "SOCK")
doParallel::registerDoParallel(workers)
message("Utilizing ",nThreads," parallel processors.")


# Parallel execution with %dopar%:
results <- foreach(i = seq(1, length(proteins))) %dopar% {
	fitLMER(fx,msstats_prot,protein=proteins[i])
}

# Close parallel connections.
suppressWarnings(stopCluster(workers))


## test contrasts -------------------------------------------------------------

# for each protein compare conditions declared in contrast_matrix
# FIXME: if run twice then errors
fit_list <- testContrasts(fit_list, contrast_matrix, moderated = TRUE)

# get results list -- df for each comparison and compute padjust
results_list <- adjustPvalues(fit_list)



