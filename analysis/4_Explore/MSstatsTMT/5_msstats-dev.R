#!/usr/bin/env Rscript

# title: SwipProteomics
# description: working through MSstats groupComparisons with repeated-measures
# design
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# data in root/rdata:
# * data_prot.rda
# * fitted.models.rda = list returned by .proposed..

## Options --------------------------------------------------------------------

## Prepare the working environment --------------------------------------------

root <- "~/projects/SwipProteomics"
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata") # temporary/large files not tracked by git


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(MSstatsTMT) # twesleyb/MSstatsTMT
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })

## load the data ---

# load msstats preprocessed protein data
load(file.path(rdatdir,"msstats_prot.rda"))
#data_prot <- msstats_prot
data <- msstats_prot

# load msstats_gene_map
load(file.path(rdatdir,"msstats_gene_map.rda"))
# gene_map

# load fitted.models
load(file.path(rdatdir,"fitted.models.rda"))

# load contrast.matrix 
load(file.path(rdatdir,"contrast.matrix.rda"))

# load dev functions
dev_fun <- list.files("./dev",pattern="*.R", full.names=TRUE)
invisible(sapply(dev_fun,source))


## do stats --------------------------------------------------------------------

# do these stats depend upon design?
# if so, then my workaround may be incorrect

# extract the linear model fitting results from list
#proteins <- fitted.models$protein # protein names
#s2.all <- fitted.models$s2 # group variance
#s2_df.all <- fitted.models$s2_df # degree freedom of s2
#lms <- fitted.models$model # linear models
#coeff.all <- fitted.models$coeff # coefficients

## init empty df to store the inference results
#num.protein <- length(proteins)
#res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) 
#colnames(res) <- c("Protein", "Comparison", "log2FC", 
#	     "pvalue", "SE", "DF", "issue")
#data$Condition <- as.factor(data$Condition) # make sure group is factor
#data$Run <- as.factor(data$Run)
#nrun <- length(unique(data$Run)) # check the number of MS runs in the data

# LOOP to perform statistical comparisons using the protein-wise models:
# ARG for testing
#i = 1
#count <- 0
#for (i in seq_along(proteins)) {
#    ## get the data for protein i

proteins = fitted.models$protein
prot = sample(proteins,1)
sub_data <- data %>% dplyr::filter(Protein == prot) 

## record the contrast matrix for each protein
sub.contrast.matrix <- contrast.matrix
sub_groups <- as.character(unique(sub_data[, c("Condition")]))

# sort the groups based on alphabetic order
sub_groups <- sort(sub_groups) 

#library(lmerTest)

## get protein's linear model
fit <- fitted.models$model[[prot]]

# actual model:
fm <- fit$model
av <- anova(fm)
s2 <- av$"Mean Sq"/av$"F value"
s2_df <- av$DenDF
coeff <- lme4::fixef(fm) 

# if not moderated s2.prior <- df.prior = 0
s2.prior = df.prior = 0
s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

# s2.post == s2 [1] TRUE # is this correct?

#s2 <- s2.all[prot] # sigma^2 # s2 <- av$"Mean Sq" / av$"F value"
#s2_df <- s2_df.all[prot] # degrees of freedom # s2_df <- av$DenDF # DF?
#coeff <- coeff.all[[prot]] # coefficients # coeff <- lme4::fixef(fit) 

#if (!is.character(fit)) { 
#s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

# LOOP to perform statistical testing for every contrast
#for (j in seq(nrow(sub.contrast.matrix))) {
#count <- count + 1
#res[count, "Protein"] <- proteins[i]
#res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] 

sub.contrast.matrix <- contrast.matrix
sub_groups <- as.character(unique(sub_data[, c("Condition")]))
sub_groups <- sort(sub_groups)

# groups with positive coefficients
j = 1; # j contrasts
idy <- sub.contrast.matrix[j, ] > 0
positive.groups <- colnames(sub.contrast.matrix)[idy]

# groups with negative coefficients
idy <- sub.contrast.matrix[j, ] < 0
negative.groups <- colnames(sub.contrast.matrix)[idy]

# make sure at least one group from each side of the contrast exist
#if (any(positive.groups %in% sub_groups) &
#  any(negative.groups %in% sub_groups)) {


.make.contrast.single(fit$model,contrast.matrix.single)


