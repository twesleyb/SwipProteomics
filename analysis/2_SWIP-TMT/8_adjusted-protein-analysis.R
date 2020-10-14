#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

stop("This script is not working when run as an executable?")

# input:
root <- "~/projects/SwipProteomics"
# * data(msstats_prot)
# * data(msstats_contrasts)

# options:
nprot = "all" # number of randomly sampled proteins to analyze
nThreads = 23 # number of cores for parallel processing
save_rda = TRUE # save results_list as rda?

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
suppressWarnings({ devtools::load_all(quiet=TRUE) })

data(gene_map)
data(partition)
data(msstats_prot)
data(msstats_contrasts)

# We need to pass the functions in lmerTestProtein to foreach, load them now:
source(file.path(root,"R","lmerTestProtein.R"))

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

# define the proteins to be analyzed
all_proteins <- unique(as.character(msstats_prot$Protein))
if (nprot == "all") {
	proteins <- all_proteins
} else {
	proteins <- sample(all_proteins, nprot)
}
message("Analyzing ",length(proteins)," protein(s).")

# loop for every protein
# Parallel execution with %dopar%
start_time <- Sys.time()
results <- list()
results <- foreach(protein = proteins, .packages=c("dplyr")) %dopar% {
	suppressPackageStartupMessages({
	  suppressMessages({
		results[[protein]] <- try(lmerTestProtein(msstats_prot,protein,fx1,cm1),silent=TRUE)
	  })
	})
}
names(results) <- proteins
stop_time <- Sys.time()

# status
message("\nElapsed time to analyze ",nprot," protein(s): ", 
	round(difftime(stop_time,start_time,units="min"),3)," (min).")

# close parallel connections
invisible(stopCluster(workers))

# remove any NULL results
drop <- which(sapply(results,is.null))
if (length(drop) > 0) {
	message("Removing ", length(drop), " NULL results.")
	results <- results[-drop]
}

# goodness of fit
# this takes a couple seconds...
r2 <- lapply(results,function(x) as.data.table(r.squaredGLMM.merMod(x$model)))
r2_df <- bind_rows(r2,.id="protein")

# Collect stats
results_df <- bind_rows(sapply(results,"[","stats")) %>% 
	left_join(r2_df, by = "protein") %>%
	filter(!isSingular) %>%
	arrange(Pvalue) %>% mutate(FDR = p.adjust(Pvalue,method="BH"))
results_df$isSingular <- NULL # drop col

# annotate with gene ids
idx <- match(results_df$protein, gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
results_df <- tibble::add_column(results_df,Symbol,.after="protein")
results_df <- tibble::add_column(results_df,Entrez,.after="Symbol")

# save the data
save(results_df, file=file.path(root,"rdata","results_df.rda"),version=2)
