#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

# input:
root <- "~/projects/SwipProteomics"
# * data(msstats_prot)
# * data(msstats_contrasts)

# options:
nprot = 1 # number of randomly sampled proteins to analyze
nThreads = 8 # number of cores for parallel processing
save_rda = FALSE # save results_list as rda?

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root,quiet=TRUE) }

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(doParallel)
})


## load the data --------------------------------------------------------------

# load msstats preprocessed protein data from SwipProteomics in root/data
devtools::load_all(quiet=TRUE)

data(msstats_prot)
data(msstats_contrasts)

# contrast matrix for Contrl-Mutant comparison:
cm1 <- setNames(c(-1,1),nm=c("GenotypeControl","GenotypeMutant"))

# lmer formula:
fx1 <- formula("Abundance ~ 0 + (1|BioFraction) + Genotype") # WT vs MUT

# the model to be fit
message("\nFitting lmer: ",fx1)

# register some nodes to do work
workers <- parallel::makeCluster(c(rep("localhost", nThreads)), type = "SOCK")
doParallel::registerDoParallel(workers)
message("\nEmploying ", nThreads," parallel processors.")

# loop for every protein
results <- list()
all_proteins <- unique(as.character(msstats_prot$Protein))

if (nprot == "all") {
	proteins <- all_proteins
} else {
	proteins <- sample(all_proteins, nprot)
}
message("Analyzing ",length(proteins)," protein(s).")

# start timer
start_time <- Sys.time()

# Parallel execution with %dopar%
results <- foreach(protein = proteins, .packages=c("dplyr","SwipProteomics")) %dopar% {
	suppressMessages({
		rho <- try(lmerTestProtein(msstats_prot,protein,fx1,cm1),silent=TRUE)
	})
	#rho$R2 <- setNames(rev(r.squaredGLMM.merMod(rho$model)),nm=c("R2c","R2fixef"))
}
stop_time <- Sys.time()

# status
message("\nElapsed time to analyze ",nprot," protein(s): ", 
	round(difftime(stop_time,start_time,units="sec"),3)," (seconds).")

if (save_rda) {
  # save the data
  message("\nSaving the data, this will take several minutes.")
  save(results,file=file.path(root,"rdata","fit1_results.rda"),version=2)
}

# close parallel connections
invisible(stopCluster(workers))
