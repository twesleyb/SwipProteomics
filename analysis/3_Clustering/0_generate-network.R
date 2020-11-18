#!/usr/bin/env Rscript 

# title: 
# description: generate protein covariation (correlation) network
# author: twab

## Input:
root <- "~/projects/SwipProteomics"

## Options:
rm_poor <- TRUE # rm proteins with overall R2 less than 0.7 before network construction

## Output:
# * generates correlation matrix which is used as input for Leidenalg clustering


## ---- prepare the working environment

# load renv
renv::load(root)

# load functions in root/R and make data in root/data accessible
# library(SwipProteomics)
devtools::load_all(root)

# load data in root/data
data(gene_map)
data(poor_prots)
data(msstats_prot)
data(msstats_results)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(neten) # twesleyb/neten
})


## ---- create covariation network

# cast the data into a matrix
dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>% as.data.table() %>% 
        as.matrix(rownames="Protein")


## examine missingness
# there are a large number or proteins with some sort of missingness
# many of these proteins are quantified in 2/3 experiments (mixtures)
# in order to retain these proteins in the network, we impute protein-
# level missing values using the KNN algorithm for MNAR data (missing-ness
# is inferred to be related to the left-shifted distribution of these 
# proteins -- they are less abundant, e.g. WASHC3)

# sum(any_missing) ~ 563 proteins with missingness
any_missing <- apply(dm,1,function(x) any(is.na(x)))

# types of missingness
n_missing <- apply(dm,1,function(x) sum(is.na(x)))

# we cant work with proteins with more than 50% missing values
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

## ---- address effect of Mixture

# we aim to identify proteins that covary together in biological space
# (BioFraction); we remove the effect of Mixture using limma prior to building
# the covariation network in order to mitigate the contribution of Mixture to
# covariation modules

# with complete cases, address effect of mixture
norm_dm <- limma::removeBatchEffect(knn_dm, batch=samples$Mixture,
			 design=model.matrix(~Condition,data=samples))

## rm poor_prots 

if (rm_poor) {
	idx <- rownames(norm_dm) %in% poor_prots
	warning("Removing ", sum(idx)," proteins before network construction.")
	norm_dm <- norm_dm[!idx,]
}

# calculate coorrelation matrix
adjm <- cor(t(norm_dm),method="pearson",use="complete.obs")


## ---- network enhancement

ne_adjm <- neten(adjm)


## ---- save networks

# coerce to data.table and save to file
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","adjm.csv")
fwrite(adjm_dt, myfile)

# coerce to data.table and save to file
ne_adjm_dt <- as.data.table(ne_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ne_adjm.csv")
fwrite(ne_adjm_dt, myfile)
