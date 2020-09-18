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
edgeR_models <- c("Common","Trended",
		  "Tagwise-Common","Tagwise-Trend")
results <- list()
for (mod in edgeR_models){
	message(paste("Simulating", mod, "dispersion model."))
	results[[mod]] <- nb.gof.m(counts=subdm, x=x, sim=n_sim,
		  model=mod, ncores=n_cores)
}

# error when using Tagwise-Common
# seemsl ike the most likely dispersion utilized by edgeR is the
# tagwise-trended dispersion

# dispersions utilized:
y <- nb.gof.m(counts=subdm, x=x, sim=n_sim,
	  model="Trended", ncores=n_cores)
summary(y)




summary(dm_gof,conv.env=0.95, data.note="Control.F4")
