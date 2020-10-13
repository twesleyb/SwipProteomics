#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: examining lmer models

# refs:
# * stargazer
#   Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary 
#   Statistics Tables. R package version 5.2.2. 
#   https://CRAN.R-project.org/package=stargazer 

## inputs 
root = "~/projects/SwipProteomics"


## prepare the environment -----------------------------------------------------
if (dir.exists(file.path(root,"renv"))) { renv::load(root) }

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(stargazer)
})

# load the required data
devtools::load_all(quiet=TRUE)
#library(SwipProteomics)
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)


## functions ------------------------------------------------------------------

subprot <- function(protein,required=c("msstats_prot")) {
	# subset the protein-level data for a given protein
	# FIXME: error message about env is not informative
	stopifnot(all(sapply(required,exists)))
	require(dplyr, quietly=TRUE)
	return(msstats_prot %>% dplyr::filter(Protein == protein))
} #EOF

submod <- function(module,required=c("msstats_prot")) {
	# subset the module-level data for a given module
	# FIXME: error message about env is not informative
	require(dplyr, quietly=TRUE)
	stopifnot(all(sapply(required,exists)))
	stopifnot("Module" %in% colnames(msstats_prot))
	return(msstats_prot %>% dplyr::filter(Module == module))
} #EOF


## main -----------------------------------------------------------------------

# annotate msstats with module membership
msstats_prot$Module <- paste0("M",partition[msstats_prot$Protein])

# get swip's module
protein <- sample(names(partition),1) 
module <- paste0("M",partition[protein])
gene <- gene_map$symbol[match(protein,gene_map$uniprot)]

# example fits
fm0 <- lmerTest::lmer("Abundance ~ (1|Mixture) + Condition",subprot(protein)) # | Condition = Genotype.BioFraction
fm1 <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + Genotype", subprot(protein))
fm2 <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + (1|Protein) + Genotype", submod(module))

# if lmerTest was used, then coerce class to something stargazer can work with
class(fm0) <- "lmerMod" 
class(fm1) <- "lmerMod" 
class(fm2) <- "lmerMod"
# NOTE: fixes error about $ and S4 class (https://stackoverflow.com/questions/31319030)

# generate latex document
cat("\\documentclass[11pt]{report}")
cat(paste0("\n\\title{",gene,"|",protein,"}"))
cat("\n\\author{TWAB}")

# title
cat("\n\\begin{document}")
cat("\n\\maketitle")

# latex tables of model summaries
cat(stargazer(fm0,title=as.character(attr(fm0,"call"))[2]))
cat(stargazer(fm1,title=as.character(attr(fm1,"call"))[2]))
cat(stargazer(fm2,title=as.character(attr(fm2,"call"))[2]))

# latex tables showing R2 of marginal and conditional types
# * marginal represents variance explained by fixed effects, e.g. Genotype
# * conditional represents the variance explained by the entire model 
#   including both (fixef and mixef).
cat(knitr::kable(r.squaredGLMM.merMod(fm0),format="latex"))
cat(knitr::kable(r.squaredGLMM.merMod(fm1),format="latex"))
cat(knitr::kable(r.squaredGLMM.merMod(fm2),format="latex"))

cat("\n\\end{document}")

# function from MuMin package
# logic from: Nakagawa, S., Schielzeth, H. (2013) A general and simple method
# for obtaining R<U+00B2> from Generalized Linear Mixed-effects Models. Methods
# in Ecology and Evolution 4: 133<U+2013>

# fixed effect is Genotype: genotype explains ~ 90% of module variance?
# mixed effects are BioFraction and Protein --> they only explain 1% of the
# variance?
