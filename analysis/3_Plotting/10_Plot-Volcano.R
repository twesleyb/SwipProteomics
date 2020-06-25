#!/usr/bin/env bash

# Plot volcano plot.
root <- getrd()
renv::load(root,quiet=TRUE)

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

devtools::load_all()

# load the data
data(tmt_protein)

plot <- ggplot()
plot <- ggplot(tmt_protein, aes(x = Adjusted.logFC, y = -log10(Adjusted.PValue)))
plot <- plot + geom_point()





