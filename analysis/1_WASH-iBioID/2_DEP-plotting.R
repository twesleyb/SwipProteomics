#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
root <- "~/projects/SwipProteomics"

# DEP se objects in bioid_se, a list that contains:
# * raw (tidy_prot) - raw bioid data 
# * sl (SL_prot) - data after sample loading norm
# * spn (SPN_prot) - data after spn norm
# * dep (dep_se) - data after imputing 
# * vsn (dep_vsn) - data after final vsn norm

## OPTIONS --------------------------------------------------------------------


## FUNCTIONS -----------------------------------------------------------------

mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


## Prepare the working environment ---------------------------------------------

# directory for output figs:
figsdir <- file.path(root,"figs","DEP"); mkdir(figsdir)


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP) # twesleyb/DEP 
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
})

# load functions in root/R
devtools::load_all(root)

# load the data in root/data
data(bioid_se) # the DEP processed data as list of se objects


## plots ----------------------------------------------------------------------
dep_se <- bioid_se$raw

p1 <- plot_frequency(dep_se)
p2 <- plot_numbers(dep_se)
p3 <- plot_coverage(dep_se)
p4 <- plot_missval(bioid_se[[1]])
p5 <- plot_detect(bioid_se[[1]])
p6 <- plot_pca(bioid_se[[1]])

p7 <- plot_normalization(bioid_se[[1]],bioid_se[[2]])
p8 <- plot_imputation(bioid_se[[1]], bioid_se[[4]])

# add final dep results (after addrejectsons) to list so we can plot:
# volcano
plot_volcano(dep, contrast = "Ubi6_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_single(dep, proteins = c("USP15", "IKBKG"))
