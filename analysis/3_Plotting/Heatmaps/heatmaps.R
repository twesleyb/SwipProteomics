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

# load the adjm (correlation matrix)
load(file.path(root,"rdata","adjm.rda"))

# network enhancment
ne_adjm <- neten::neten(adjm,alpha=0.9,diffusion=0.4)

# all modules
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# examine a random module
m <- sample(names(modules),1)
m = "M4"
prots <- modules[[m]]
sub_adjm <- adjm[prots,prots]
sub_ne_adjm <- ne_adjm[prots,prots]
colnames(sub_adjm) <- colnames(sub_ne_adjm) <- mapID(prots,"uniprot","symbol")
p0 <- ggcorrplot(sub_adjm, hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c(wt_color,"gray",mut_color),
		   type = "upper")
# generate plots
p1 <- ggcorrplot(sub_ne_adjm, hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c("gray","gray",mut_color),
		   type = "lower")
p0 <- p0 + ggtitle(m)
cowplot::plot_grid(p0,p1)

# save figures
myfile <- file.path(figsdir,"adjm.pdf")
ggsave(myfile,p0,width=5,height=5)

myfile <- file.path(figsdir,"ne_adjm.pdf")
ggsave(myfile,p1,width=5,height=5)
