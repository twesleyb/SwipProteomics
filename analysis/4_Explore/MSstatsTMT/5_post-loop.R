#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: collect bash loop results for WTvMutant comparisons
# combine with intrafraction comparisons and save

# input:
# * output from loop, saved as data/loop_results.rda
# * msstats_gene_map

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

devtools::load_all()

# load gene map
myfile <- file.path(root,"rdata","msstats_gene_map.rda")
load(myfile)

# load results from loop 
data(loop_results)
# results

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
message("Summary of significant proteins (Control-Mutant):")
df <- data.table::data.table("N sig (FDR<0.05)" = sum(results$FDR < 0.05),
		 "N sig (PAdj<0.05)" = sum(results$PAdjust < 0.05))
knitr::kable(df)

# top sig proteins
message("Top 5 significant proteins:")
results %>% arrange(PValue) %>% select(Symbol,log2FC,PValue,PAdjust,Tstatistic,SE,DF) %>% head(5) %>% knitr::kable()


# combine with other results results -------------------------------------------

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

# drop NA
results_list <- lapply(results_list,function(x) x[!is.na(x$PValue),])

# summary
message("Summary of significant proteins (Control.F#-Mutant.F#):")
knitr::kable(t(sapply(results_list,function(x) sum(x$PAdjust < 0.05))))

# list of sig prots
sig_prots <- lapply(results_list, function(x) {
			    x %>% filter(PAdjust < 0.05) %>% 
				    select(Protein) %>% 
				    unlist() %>% 
				    as.character() })

all_prots <- Reduce(union,sig_prots)
message("Total number of unique sig prots: ", length(all_prots))

# total
overlap <- Reduce(intersect,sig_prots)
names(overlap) <- gene_map$symbol[match(overlap,gene_map$uniprot)]

# overlap
message("Commonly significant sig prots: ")
knitr::kable(overlap)

# save
final_results = c(results_list[c("F4","F5","F6","F7","F8","F9","F10")],"Control-Mutant"=list(results))
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(final_results,myfile)


## Compare this to our previous results ----------------------------------------

# load the swip data
data(swip_tmt)

# edgeR intrafraction
df0 <- swip_tmt %>% select(Accession,Fraction,PValue,FDR,logFC) %>% unique()
colnames(df0) <- c("protein","contrast","pval","fdr","log2fc")

# edgeR summary
message("Summary of edgeR sig prots (FDR<0.05) for intrafraction contrasts:")
df0 %>% group_by(contrast) %>% summarize(n=sum(fdr<0.05)) %>% knitr::kable()

# edgeR WTvMut
df1 <- swip_tmt %>% select(Accession,Adjusted.PValue,Adjusted.FDR,Adjusted.logFC) %>% unique()
colnames(df1) <- c("protein","pval","fdr","log2fc")
df1$contrast <- "Mutant-Control"

# edge summary
message("Summary of EdgeR sig prots for Control-Mutant contrast:")
df1 %>% group_by(contrast) %>% 
	summarize(`FDR < 0.05`=sum(fdr<0.05)) %>% knitr::kable()

# MSstats intrafraction
df2 <- msstats_results %>% filter(!is.na(PValue)) %>%
       	select(Protein, Contrast, PValue, log2FC)
colnames(df2) <- c("protein","contrast","pval","log2fc")
df2$contrast <- gsub("Mutant.F[0-9]{1,2}-Control.","",df2$contrast)

# edgeR WTvMut
df3 <- results %>% filter(!isSingular) %>% 
	filter(!is.na(PValue)) %>% 
	select(Protein, Contrast, PValue,log2FC)
colnames(df3) <- c("protein","contrast","pval","log2fc")

# combine intrafraction restults
df_x = left_join(df0,df2,by=c("protein","contrast"),suffix=c(".edgeR",".MSstats"))

# combine WTvMut results
df_y = left_join(df1,df3,by=c("protein","contrast"),suffix=c(".edgeR",".MSstats"))

# coorelation of pvals
rho0 <- cor(df_x$pval.edgeR,df_x$pval.MSstats,method="spearman",use="pairwise.complete")
rho1 <- cor(df_y$pval.edgeR,df_y$pval.MSstats,method="spearman",use="pairwise.complete")

message("Rank-Correlation of PValues:")
data.table("Control.F#-Mutant.F#" = rho0, "Control-Mutant" = rho1) %>%
	knitr::kable()

# coorelations of fold change
rho0 <- cor(df_x$log2fc.edgeR,df_x$log2fc.MSstats,method="spearman",use="pairwise.complete")
rho1 <- cor(df_y$log2fc.edgeR,df_y$log2fc.MSstats,method="spearman",use="pairwise.complete")

message("Rank-Correlation of log2fc:")
data.table("Control.F#-Mutant.F#" = rho0, "Control-Mutant" = rho1) %>%
	knitr::kable()

# the results are highly correlated
