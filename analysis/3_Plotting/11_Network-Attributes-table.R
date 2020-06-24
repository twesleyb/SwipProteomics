#!/usr/bin/env Rscript

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports.
devtools::load_all()
suppressPackageStartupMessages({
	library(ggplot2)
	library(grid)
	library(gridExtra)
})

# Set plotting theme.
set_font("Arial",font_path=file.path(root,"fonts"))

# Load the data.
data(partition)
data(module_noa)
data(tmt_protein)

# Summarize key module statistics.
df0 <- data.frame("Nodes" = formatC(length(unique(tmt_protein$Accession)),big.mark=","),
		 "Modules" = length(unique(partition))-1,
		 "Percent Clustered" = round(100*sum(partition!=0)/length(partition),2))
df1 <- module_noa %>% filter(Module != "M0") %>% 
	summarize("Median Size" = median(`Module Size`),
		  "Median PVE" = round(100*median(`Module PC1 PVE`),2))

# Create table.
gtab <- tableGrob(cbind(df0,df1),rows=NULL, theme=ttheme_default())

# Save.
myfile <- file.path(root,"figs","Networks","Network_Summary.pdf")
ggsaveTable(gtab,myfile)
