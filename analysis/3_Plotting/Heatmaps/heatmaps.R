#!/usr/bin/env Rscript

## inputs 
root = "~/projects/SwipProteomics"
wt_color = "#47b2a4"

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
if (!dir.exists(figsdir)) { dir.create(figsdir); message("mkdir ",figsdir) }

# load the required data
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# load data in root rdata -- to big to be stored in root/data
load(file.path(root,"rdata","adjm.rda"))
#load(file.path(root,"rdata","ne_adjm.rda"))
#load(file.path(root,"rdata","ppi_adjm.rda"))

# collect list of modules
modules = split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

#for (module in modules){
module = "M25"
prots <- modules[[module]]
prots <- prots[prots %in% colnames(adjm)]


p0 <- ggcorrplot(ne_adjm[prots,prots], 
		   hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c(module_colors[[module]],"gray",wt_color),
		   type = "lower")

p1 <- ggcorrplot(adjm[prots,prots], 
		   hc.order = TRUE, 
		   outline.col = "white", 
		   colors = c(module_colors[[module]],"gray",wt_color),
		   type = "upper")

p0 <- p0 + ggtitle(paste(module,"mean(bicor) =",mean(subadjm0)))

# Customize x-axis label color
build <- ggplot_build(p0)
x_labels <- build$layout$panel_params[[1]]$x$get_labels()
sig_prots <- msstats_results %>% filter(Protein %in% prots) %>%
	filter(adj.pvalue < 0.05) %>% select(Protein) %>% unique() %>% unlist()
x_colors <- rep("black",ncol(subadjm0))
x_colors[x_labels %in% sig_prots] <- "dark red"
p0 <- p0 + theme(axis.text.x = element_text(color = x_colors))

# need to specify type both and then  edit data?
