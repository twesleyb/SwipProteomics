#!/usr/bin/env Rscript

# title: SwipProteomics
# description: anaysis of the goodness of fit of negative binomial models
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Misc functions -------------------------------------------------------------

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

# Prepare the R environment ---------------------------------------------------

# load renv
root <- getrd()
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(NBGOF)
	library(dplyr)
	library(data.table)
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# load data from NBGOF package
data(swip_tmt) # tmt_protein


# consider the data from a single group, 3x biological replicates
dm <- swip_tmt %>%
	filter(Treatment == "Control", Fraction == "F4") %>%
	dcast(Accession ~ Sample, value.var = "Intensity") %>%
	as.data.table() %>%
	as.matrix(rownames="Accession")
idx <- sample(nrow(dm),500) # sample 500 rows
subdm <- dm[idx,]

# number of cores for parallel processing
n_cores <- parallel::detectCores() - 1

# simple model matrix for single group, three replicates
x <- as.matrix(rep(1,3))

# GOF tests for different dispersion models:
# edgeR models: Common, Trended, Tagwise-Common,Tagwise-Trend
n_sim = 999
models <- c("Common"="CoxReid",
	    "Trended"="auto",
	    "Genewise"="auto")
results <- list()
for (mod in names(models)){
	message(paste("Simulating", mod,
		      "dispersion model."))
	gof <- nb.gof.m(counts=subdm,
			x=as.matrix(rep(1,3)),
			sim=n_sim,
		        model=mod,
		        method=models[mod],
			ncores=n_cores)
	summary(gof,conv.env=0.95,data.node="SWIP TMT")
	results[[mod]] <- gof
}

# save the restults
myfile <- file.path(root,"rdata","dispersion_simulations.RData")
saveRDS(results,file=myfile)
