#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: 

# input:
# * msstats_prot.rda

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
})

# load gene map
myfile <- file.path(root,"rdata","msstats_gene_map.rda")
load(myfile)

# load statistical results
results_df <- data.table::fread("results.csv")

# Padjust
padj <- p.adjust(results_df$Pvalue,method="BH")
results_df <- tibble::add_column(results_df,"Padjust"=padj,.after="Pvalue")

# annotate results with gene IDs
accession <- results_df$protein
idx <- match(accession,gene_map$uniprot)
symbol <- gene_map$symbol[idx]
entrez <- gene_map$entrez[idx]
results_df <- tibble::add_column(results_df, symbol, .after="protein")
results_df <- tibble::add_column(results_df, entrez, .after="symbol")

results_df %>% filter(Padjust < 0.05) %>% arrange(Pvalue) %>% head(50) %>% knitr::kable()
