#!/usr/bin/env Rscript

## inputs 
root = "~/projects/SwipProteomics"

# explore ne params

## functions ------------------------------------------------------------------


## main ------------------------------------------------------------------------

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root) }

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(ggcorrplot)
})

# load project
devtools::load_all(quiet=TRUE)

# dir for output
figsdir <- file.path(root,"figs","Heatmaps")
if (!dir.exists(figsdir)) { dir.create(figsdir) }

# load the required data
data(swip)
data(gene_map)
data(wt_color)
data(mut_color)
data(partition)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# load data in root rdata -- to big to be stored in root/data
load(file.path(root,"rdata","adjm.rda"))
load(file.path(root,"rdata","ne_adjm.rda"))
load(file.path(root,"rdata","ppi_adjm.rda"))

# prepare the data
x = adjm[prots,prots]
y = ne_adjm[prots,prots]
z = y*(1/max(y)) # scale such that max is 1

## explore neten params
# alpha = small ==> super sparse
# diffusion = 2 == most sparse
# use alpha to adjust the overall connectivity -- smaller then sparser
y = neten::neten(adjm,alpha=0.9,diffusion=1.0) 
x = y[prots,prots]
p0 <- ggcorrplot(x, hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c("gray","gray",mut_color),
		   type = "upper")
p0

colnames(x) <- colnames(z) <- mapID(prots,"uniprot","symbol")

p0 <- ggcorrplot(x, hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c(wt_color,"gray",mut_color),
		   type = "upper")

p1 <- ggcorrplot(z, hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c("gray","gray",mut_color),
		   type = "lower")


# save figures
myfile <- file.path(figsdir,"adjm.pdf")
ggsave(myfile,p0,width=5,height=5)

myfile <- file.path(figsdir,"ne_adjm.pdf")
ggsave(myfile,p1,width=5,height=5)
