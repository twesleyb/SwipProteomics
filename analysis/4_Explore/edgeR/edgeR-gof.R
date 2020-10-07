#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: assess the goodness-of-fit of SWIP proteomics data with edgeR
#' authors: Tyler W A Bradshaw
#' ---

## INPUT:
# * swip_tmt the preprocessed SWIP proteomics dataset

## OPTIONS:

## OUTPUT:
# * plots in root/figs/edgeR

## Prepare the R working environment ------------------------------------------

root <- "~/projects/SwipProteomics/"
renv::load(root,quiet=TRUE)

# load required packages
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(data.table) # for working with tables
	library(limma) # required by edgeR
	library(edgeR) # for statistical analysis
	#library(fitdistrplus) # for fitting nb
})

# load project specific functions and data
suppressWarnings({ devtools::load_all() })

# load the data
data(samples) # swip tmt sample meta data
data(swip_tmt) # swip tmt preprocessed protein data


## pass data to edgeR ---------------------------------------------------------

# cast the normalized data into a matrix
dm <- swip_tmt %>% filter(Treatment != "SPQC") %>%
	dcast(Accession ~ Sample, value.var="Intensity") %>%
	as.matrix(rownames="Accession")

# create dge object
dge <- edgeR::DGEList(counts=dm)

# perform library (run-level) normalization
dge <- calcNormFactors(dge,method="TMM")

# names(dge)
#[1] "counts" "samples"

# map sample meta data onto dge object
sample_names <- rownames(dge$samples)
idx <- match(sample_names,samples$Sample)
dge$samples$Fraction <- samples$Fraction[idx]
dge$samples$Genotype <- samples$Genotype[idx]
# update sample group(s) -- Fraction.Genotype
dge$samples$group <- interaction(samples$Fraction[idx],samples$Genotype[idx])

# create a design matrix with update samples
#design <- model.matrix(~ 0 + group, data=dge$samples)
design <- model.matrix(~ Fraction + Genotype, data=dge$samples)


## estimate dispersion ---------------------------------------------------------
# edgeR supports several methods of estimating dispersion.
# For fitting GLMs it is only appropriate to use a global metric such as
# 'Common' or 'Trended' dispersion. Here we explicitly estimate 'Trended' 
# dispersion to account for Mean - Variance relationship in protein
# quantification.
#dge <- estimateGLMTrendedDisp(dge, design, method="auto")

# estimate all dispersion using all methods 
dge <- estimateDisp(dge,design,method="auto")

#dge <- estimateGLMCommonDisp(dge,design,method="CoxReid")

# NOTE: edgeR also implements a glmFit and glmLRT test. The glmQLFit and
# glmQLFTest are more robust. Only common, trended dispresion can be used.

## fit protein-wise NB GLMs ---------------------------------------------------
# fit model and generate data QC plots

fit <- glmQLFit(dge, design, robust = TRUE)

# BCV
myfile <- file.path(root,"figs","edgeR","plotBVC.pdf")
pdf(file=myfile)
plotBCV(dge)

#clean-up
invisible(dev.off())

# QLDispersion
myfile <- file.path(root,"figs","edgeR","plotQLDisp.pdf")
pdf(file=myfile)
plotQLDisp(fit)

#clean-up
invisible(dev.off())


## assess GOF -----------------------------------------------------------------
# use edgeR gof function
# gof
myfile <- file.path(root,"figs","edgeR","gof.pdf")
pdf(file=myfile)
edgeR_gof <- edgeR::gof(fit,plot=TRUE)

# clean-up
invisible(dev.off())

message("\nSummary of outliers (blue):")
knitr::kable(table(edgeR_gof$outlier))
idx <- edgeR_gof$outlier
to_drop  <- rownames(dge)[idx]

# DONE
