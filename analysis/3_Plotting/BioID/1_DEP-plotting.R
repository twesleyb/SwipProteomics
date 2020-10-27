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

FDR_alpha = 0.1
enrichment_threshold = log2(3.0)

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
        library(ggplot2) # for manipulating plots
	library(data.table) # for working with data.tables
})

# load functions in root/R
devtools::load_all(root)

# load the data in root/data
data(bioid_se) # the DEP processed data as list of se objects

# names(bioid_se):
# [1] "raw"     "sl"      "spn"     "dep"     "vsn"     "results"

## plots ----------------------------------------------------------------------

raw <- bioid_se$raw # raw data
sln <- bioid_se$sl # post-SL norm
spn <- bioid_se$spn # post-SPN norm
imp <- bioid_se$dep # post-imputing
vsn <- bioid_se$vsn # post-vsn
dep <- bioid_se$results

#dep <- add_rejections(dep, alpha = FDR_alpha, lfc = enrichment_threshold)

df <- DEP::get_results(dep)

yint <- -1*log10(max(df$WASH_vs_Control_p.val[df$WASH_vs_Control_p.adj<FDR_alpha]))

myfile <- file.path(figsdir,"WASH_BioID_DEP_plots.pdf")
pdf(file=myfile,onefile=TRUE)

  hist(df$WASH_vs_Control_p.val,breaks=50)
  abline(v=y)

  #meanSdPlot(raw)

  #meanSdPlot(vsn)

  #plot_frequency(raw)

  #plot_numbers(raw)

  #plot_coverage(raw)

  #plot_missval(raw)

  plot_detect(raw)

  plot_normalization(raw,vsn)

  plot_pca(vsn)

  #plot_imputation(spn, imp)

  plot <- plot_volcano(dep, contrast = "WASH_vs_Control", 
		     label_size = 2, add_names = TRUE)
  plot <- plot + geom_vline(xintercept=enrichment_threshold,linetype="dashed",color="darkred")
  plot <- plot + geom_hline(yintercept=yint,linetype="dashed",color="darkred")
  plot


invisible(dev.off())
