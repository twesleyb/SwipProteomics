#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description:
#' authors: Tyler W A Bradshaw
#' ---

## INPUT:
# data(swip_tmt) # the preprocessed SWIP proteomics dataset

## OPTIONS:

## OUTPUT:

## FUNCTIONS ------------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# Get the repository's root directory.
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


## MAIN prepare the R working environment -------------------------------------

root <- getrd()
renv::load(root,quiet=TRUE)

# load required packages
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(data.table) # for working with tables
	library(limma) # required by edgeR
	library(edgeR) # for statistical analysis
})

# load project specific functions and data
suppressWarnings({ devtools::load_all() })

# load the data
data(samples) # swip sample meta data
data(swip_tmt) # swip tmt


## ----------------------------------------------------------------------------

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
design <- model.matrix(~ 0 + group, data=dge$samples)

# estimate dispersion
# edgeR supports several methods of estimating dispersion
# for GLM it is only appropriate to use Common or Trended
# dispersion. Explicitly estimate Trended dispersion to
# account for Mean - Variance relationship.
#dge <- estimateDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design, method="auto")

# NOTE: edgeR also implements a glmFit and glmLRT test. The glmQLFit and
# glmQLFTest are more robust. Only common, trended dispresion can be used.

fit <- glmQLFit(dge, design, robust = TRUE)
