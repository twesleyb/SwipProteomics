#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of intra-fraction differential abundance with DEP
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify projject's root directory
ROOT <- "~/projects/SwipProteomics"

## OPTIONS --------------------------------------------------------------------
FDR_alpha <- 0.05 # BH significance threshold for protein DA

## FUNCTIONS -----------------------------------------------------------------

fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}


mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


## Prepare the working environment ----------------------------------------------

# directory for output tables:
tabsdir <- file.path(ROOT,"tables"); mkdir(tabsdir)
datadir <- file.path(ROOT,"data"); mkdir(datadir)


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(ROOT,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP) # twesleyb/DEP
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
})

# load functions in root/R
devtools::load_all(ROOT)

# load the data in root/data
data() # raw data
data(swip_tmt)
data(samples)
data(gene_map)


## fix gene_map ----------------------------------------------------------------

# fix missing gene name in gene_map
idx <- which(is.na(gene_map$symbol))
if (gene_map$uniprot[idx] == "Q80WG5") {
	gene_map$symbol[idx] <- "Lrrc8a"
}


