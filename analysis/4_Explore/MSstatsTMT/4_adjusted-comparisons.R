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


## Perform protein-level statistical tests --------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment.

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


#------------------------------------------------------------------------------
## do stats

# do these stats depend upon design?
# if so, then my workaround may be incorrect

# extract the linear model fitting results from list
proteins <- fitted.models$protein # protein names
s2.all <- fitted.models$s2 # group variance
s2_df.all <- fitted.models$s2_df # degree freedom of s2
lms <- fitted.models$model # linear models
coeff.all <- fitted.models$coeff # coefficients

## init empty df to store the inference results
num.protein <- length(proteins)
#res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) 
#colnames(res) <- c("Protein", "Comparison", "log2FC", 
#	     "pvalue", "SE", "DF", "issue")
data$Condition <- as.factor(data$Condition) # make sure group is factor
data$Run <- as.factor(data$Run)
nrun <- length(unique(data$Run)) # check the number of MS runs in the data

# LOOP to perform statistical comparisons using the protein-wise models:
# ARG for testing
#i = 1
#count <- 0
#for (i in seq_along(proteins)) {
#    ## get the data for protein i

i = 1
prot = proteins[i]
sub_data <- data %>% dplyr::filter(Protein == prot) 

## record the contrast matrix for each protein
sub.contrast.matrix <- contrast.matrix
sub_groups <- as.character(unique(sub_data[, c("Condition")]))

# sort the groups based on alphabetic order
sub_groups <- sort(sub_groups) 

#library(lmerTest)

## get protein's linear model
model_list <- lms[[prot]] 

# actual model:
fm <- model_list$model
av <- anova(fm)
s2 <- av$"Mean Sq"/av$"F value"
s2_df <- av$DenDF
coeff <- lme4::fixef(fm) 
# if not moderated, s2.prior <- df.prior = 0
s2.prior = df.prior = 0
s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
# s2.post == s2 [1] TRUE # is this correct?

s2 <- s2.all[prot] # sigma^2 # s2 <- av$"Mean Sq" / av$"F value"
s2_df <- s2_df.all[prot] # degrees of freedom # s2_df <- av$DenDF # DF?
coeff <- coeff.all[[prot]] # coefficients # coeff <- lme4::fixef(fit) 

#if (!is.character(fit)) { 
s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

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

###############################################################################

contrast.matrix.single <- as.vector(sub.contrast.matrix[j, ])
names(contrast.matrix.single) <- colnames(sub.contrast.matrix)
#cm = contrast1
cm <- setNames(c(0,1),nm=c("(Intercept)", "ConditionMutant"))
FC <- (cm %*% coeff)[, 1]

## variance and df
#if (inherits(fm, "lm")) {
# lm
#variance <- diag(t(cm) %*% 
#	     summary(fm)$cov.unscaled %*% cm) * s2.post
#df.post <- s2_df + df.prior

# lmer
# returns a function that takes Lc and thpars
vss <- .vcovLThetaL(fm)
isLMM(fm) # TRUE
#calling update() to refit model

vss <- function(fm, Lc, thpars) {
	#fm2 = update(fm,REML=FALSE)# do we need to set REML?
	#stopifnot(is.numeric(thpars), length(thpars) == np)
	np = length(fm@pp$theta)
	sigma2 <- thpars[np]^2
	#ff2 = update(fm, devFunOnly = TRUE) # @ function theta
	ff2(thpars[-np])
	ff2 is function?
	# what is RXi? pp is of class merPredD from lme4- some sort of generator
	tcrossprod(attr(fm,"pp")$RXi()) # access RXi within fm
	vcov_unscaled <- tcrossprod(envff2$pp$RXi())
	vcov_out <- sigma2 * vcov_unscaled
	return(list(
	varcor = as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)),
	unscaled.varcor = vcov_unscaled,
	sigma2 = sigma2
} #class(ans) <- ".vcovLThetaL"


## for the theta and sigma parameters:
# (Intercept) ConditionMutant
varcor <- vss(t(cm), c(fit$thopt, fit$sigma)) 

            vcov <- varcor$unscaled.varcor * s2
            se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)

            ## calculate variance
            vcov.post <- varcor$unscaled.varcor * s2.post
            variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)

            ## calculate df see [14]
            g <- .mygrad(function(x) vss(t(cm),x)$varcor,c(fit$thopt,fit$sigma))
            denom <- try(t(g) %*% fit$A %*% g, silent = TRUE)
            if (inherits(denom, "try-error")) {
              df.post <- s2_df + df.prior
            } else {
              df.post <- 2 * (se2)^2 / denom + df.prior
            }
          } # EIS
          ## calculate the t statistic
          t <- FC / sqrt(variance)
          ## calculate p-value
          p <- 2 * pt(-abs(t), df = df.post)
          res[count, "pvalue"] <- p
          ## save testing results
          res[count, "log2FC"] <- FC
          res[count, "SE"] <- sqrt(variance)
          res[count, "DF"] <- df.post
          if (s2_df == 0) {
            res[count, "issue"] <- "SingleMeasurePerCondition"
          } else {
            res[count, "issue"] <- NA
          }
          # continued from:
	  # if (any(positive.groups %in% sub_groups) &
          #  any(negative.groups %in% sub_groups)) {
        } else {
          # at least one condition is missing
          out <- .issue.checking(
            data = sub_data,
            contrast.matrix = sub.contrast.matrix[j, ]
          )
          res[count, "log2FC"] <- out$logFC
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          res[count, "issue"] <- out$issue
        } #EIS
      } # ENDS LOOP
