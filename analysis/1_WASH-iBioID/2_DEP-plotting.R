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

# names(bioid_se):
# [1] "raw"     "sl"      "spn"     "dep"     "vsn"     "results"

## plots ----------------------------------------------------------------------

raw <- bioid_se$raw
sln <- bioid_se$sl
spn <- bioid_se$spn
imp <- bioid_se$dep
vsn <- bioid_se$vsn
dep <- bioid_se$results

# generate plots
suppressMessages({
  suppressWarnings({

myfile <- file.path(figsdir,"WASH_BioID_DEP_plots.pdf")
pdf(file=myfile,onefile=TRUE)
  print(meanSdPlot(vsn))
  print(plot_frequency(raw))
  print(plot_numbers(raw))
  print(plot_coverage(raw))
  print(plot_missval(raw))
  print(plot_detect(raw))
  print(plot_normalization(raw,vsn))
  print(plot_pca(vsn))
  print(plot_imputation(spn, imp))
  print(plot_volcano(dep, contrast = "WASH_vs_Control", label_size = 2, add_names = TRUE))
invisible(dev.off())

  })
})

# NOTE: boxplot function seems be be broken
