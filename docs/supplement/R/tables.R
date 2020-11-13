#!/usr/bin/env Rscript

# title: supplment.R
# description: code associated with supplement.Rnw
# author: twab

# NOTE: code chunks are delimited with '## ---- LABEL' where LABEL cooresponds
# to its label in the main *.Rnw file. Declare code chunk options in the main
# R noweave file. 

# NOTE: it seems best if your approach each chunk as a stand-alone bit of code.

# prepare the renv
root <- "~/projects/swipproteomics"
renv::load(root)
devtools::load_all(root)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- tab0

library(grid)
library(gtable)
library(gridExtra)
library(ggplot2)
library(cowplot)

df <- summary(fm0, ddf="Satterthwaite")[["coefficients"]] %>%
	as.data.table(keep.rownames="Coefficient")
df$Coefficient <- gsub("Genotype|BioFraction","",df$Coefficient)
colnames(df)[colnames(df)=="Pr(>|t|)"] <- "p value"
colnames(df)[colnames(df)=="Std. Error"] <- "SE"
colnames(df)[colnames(df)=="df"] <- "DF"
df$"p value" <- formatC(df$"p value",digits=3)
df$Estimate <- round(df$Estimate,2)
df$SE <- round(df$SE,3)
df$DF <- round(df$DF,2)
df$"t value" <- round(df$"t value",2)
tt <- gridExtra::ttheme_default(
      base_size = 11,
      core = list(bg_params = list(fill = "white"))
    )
tab <- tableGrob(df, rows = NULL, theme = tt)
# add outside border
g <- gtable_add_grob(tab, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), 
		     t = 1, b = nrow(tab), l = 1, r = ncol(tab))
# add rows -- i have no idea how gtable works
for (i in c(1:nrow(tab))) {
	g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), 
		     t = i, b = i, l = 1, r = ncol(tab))
}
plot_grid(g)
