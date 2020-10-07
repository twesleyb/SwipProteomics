#!/usr/bin/env Rscript

# title: SwipProteomics
# description: using variancePartition and mixed effect linear models for
#   hypothesis testing
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
REPO <- "SwipProteomics"
ROOT <- file.path("~/projects/", REPO)

## OPTIONS --------------------------------------------------------------------
ncores <- -1

## Prepare the R environment --------------------------------------------------

# load renv
renv::load(ROOT,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(variancePartition) # for mixed effects linear modeling
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
	library(BiocParallel) # for parallel processing
})

# load projects's data and functions
suppressPackageStartupMessages({ devtools::load_all(ROOT) })

# load the data in root/data
data(swip_tmt)
data(samples)
data(gene_map)

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
if (ncores == -1) { 
	n = detectCores()-1 
} else {
	n = ncores
}
param = SnowParam(n, "SOCK", progressbar=TRUE)
register(param)

## the experimental design ----------------------------------------------------

# insure levels of Fraction are in order
levels(samples$Fraction) <- c("F4","F5","F6","F7","F8","F9","F10")

# the experimental design:
exp_design <- samples %>% select(Sample,Experiment,Channel,Treatment,Fraction) %>% 
	arrange(Experiment,Treatment)
rownames(exp_design) <- exp_design$Sample
knitr::kable(exp_design,row.names=FALSE)


## fix gene_map ---------------------------------------------------------------

# fix missing gene name in gene_map
idx <- which(is.na(gene_map$symbol))
if (gene_map$uniprot[idx] == "Q80WG5") {
	gene_map$symbol[idx] <- "Lrrc8a"
}


## Mixed Effects Linear Models ------------------------------------------------

# The variable to be tested must be a fixed effect.
# NOTE, this contrasts from the model from above: 
# model.matrix( ~ Disease, metadata)
# We now pass the blocking factor Individual as a mixed effect.
form <- ~ 0 + Treatment + (1|Fraction)

# estimate weights using linear mixed model of dream
voom = voomWithDreamWeights(prot_dm, form, exp_design)

L = getContrast(voom, form, exp_design, c("TreatmentControl", "TreatmentMutant"))

print(L) # NOTE: Would a matrix be produced if design was more complex?

# Visualize contrast matrix
#plotContrasts(L) 

# fit dream model with contrasts
fit = dream(voom, form, exp_design, L)

# get names of available coefficients and contrasts for testing
colnames(fit)

# extract results from first contrast
topTable(fit, coef="TreatmentMutant", number=3)
