#!/usr/bin/env Rscript 

# title: 
# description: generate protein covariation (correlation) network
# author: twab

## Input:
root <- "~/projects/SwipProteomics"


## Options

## prepare the working environment ---------------------------------------------

# load renv
renv::load(root)

# load functions in root/R and make data in root/data accessible
# library(SwipProteomics)
devtools::load_all(root)

# load data in root/data
data(gene_map)
data(msstats_prot)
data(msstats_results)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## create covariation network -------------------------------------------------

# cast the data into a matrix
dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>% as.data.table() %>% 
        as.matrix(rownames="Protein")


# 563 proteins with missingness
any_missing <- apply(dm,1,function(x) any(is.na(x)))

# flag proteins with missingness--we will not use these in the statistical
# analysis of modules.

## types of missingness
n_missing <- apply(dm,1,function(x) sum(is.na(x)))

# drop proteins with more than 50% missingness
drop <- apply(dm,1,function(x) sum(is.na(x))> 0.5 * ncol(dm))

# knn impute
knn_data <- impute::impute.knn(dm[!drop,])
knn_dm <- knn_data$data

# build sample metadata
namen <- colnames(knn_dm)
samples <- as.data.table(do.call(rbind,strsplit(namen,"_")))
colnames(samples) <- c("Mixture","Genotype","BioFraction")
samples$Condition <- interaction(samples$Genotype,samples$BioFraction)

# with complete cases, address effect of mixture
norm_dm <- limma::removeBatchEffect(knn_dm, batch=samples$Mixture,
			 design=model.matrix(~Condition,data=samples))

# calculate coorrelation matrix
adjm <- cor(t(norm_dm),method="pearson",use="complete.obs")

# coerce to data.table and save to file
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","adjm.csv")
fwrite(adjm_dt, myfile)
