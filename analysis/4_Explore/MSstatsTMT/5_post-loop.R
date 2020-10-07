#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: 

# input:

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
})

devtools::load_all()

# load gene map
myfile <- file.path(root,"rdata","msstats_gene_map.rda")
load(myfile)

# load results from loop 
results <- data.table::fread("results.csv")

# compute FDR and Padjust
fdr <- p.adjust(results$Pvalue,method="BH")
padj <- p.adjust(results$Pvalue,method="bonferroni")
results <- tibble::add_column(results,"FDR"=fdr,.after="Pvalue")
results <- tibble::add_column(results,"Padjust"=padj,.after="FDR")

# annotate results with gene symbols
accession <- results$protein
idx <- match(accession,gene_map$uniprot)
symbol <- gene_map$symbol[idx]
entrez <- gene_map$entrez[idx]
data.table::setnames(results, "protein", "Protein")
results <- tibble::add_column(results, Symbol=symbol, .after="Protein")
results <- tibble::add_column(results, Entrez=entrez, .after="Symbol")

# Munge names  - run once!
data.table::setnames(results, "Padjust", "PAdjust")
data.table::setnames(results, "Pvalue", "PValue")
data.table::setnames(results, "contrast", "Contrast")

# summarize sig prots
df <- data.table::data.table("N sig (FDR<0.05)" = sum(results$FDR < 0.05),
		 "N sig (PAdj<0.05)" = sum(results$PAdjust < 0.05))
knitr::kable(df)

# top sig proteins
results %>% arrange(PValue) %>% select(Symbol,log2FC,PValue,PAdjust,Tstatistic,SE,DF) %>% head(5) %>% knitr::kable()


# combine with intrafractino results -------------------------------------------

# load the data
load(file.path(root,"rdata","msstats_results.rda"))

# annotate results with gene symbols
accession <- msstats_results$Protein
idx <- match(accession,gene_map$uniprot)
symbol <- gene_map$symbol[idx]
entrez <- gene_map$entrez[idx]
msstats_results <- tibble::add_column(msstats_results, Symbol=symbol, .after="Protein")
msstats_results <- tibble::add_column(msstats_results, Entrez=entrez, .after="Symbol")

# munge harder
contrasts <- msstats_results$Label
msstats_results <- tibble::add_column(msstats_results, Contrast=contrasts, .before="Protein")
msstats_results$Label <- NULL

# Munge names - run once!
data.table::setnames(msstats_results, "pvalue", "PValue")
data.table::setnames(msstats_results, "adj.pvalue", "PAdjust")

# split msstats results into a list by Contrast (or Fraction)
results_list <- msstats_results %>% group_by(Contrast) %>% group_split()
names(results_list) <- sapply(results_list,function(x) unique(x$Contrast))
names(results_list) <- gsub("Mutant.F[0-9]{1,2}-Control.","",names(results_list))

# save
final_results <- list("Control-Mutant"=list(results),results_list[c("F4","F5","F6","F7","F8","F9","F10")])
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(final_results,myfile)
