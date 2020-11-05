#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(swip)
data(swip_tmt) # edgeR results
data(msstats_results)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

edgeR_df <- swip_tmt %>% select(Accession,Adjusted.logFC,Adjusted.PValue) %>% 
	unique()
colnames(edgeR_df) <- c("Protein","log2FC","Pvalue")
edgeR_df$Method <- "edgeR"

MSstatsTMT_df <- msstats_results %>% select(Protein,log2FC,Pvalue)
MSstatsTMT_df$Method <- "MSstatsTMT"

df <- full_join(edgeR_df, MSstatsTMT_df, 
		by=c("Protein","Pvalue","log2FC","Method"))
df$Method <- factor(df$Method, levels=c("edgeR","MSstatsTMT"))

library(ggplot2)

plot <- ggplot(data=df)
plot <- plot + aes(x=Pvalue, fill=Method, group=Method)
plot <- plot + geom_histogram(bins=50)
