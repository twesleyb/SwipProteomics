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
  library(ggplot2)
})

devtools::load_all()

# load gene map
data(msstats_gene_map)

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
results %>% 
	arrange(PValue) %>% 
	select(Symbol,log2FC,PValue,PAdjust,Tstatistic,SE,DF) %>% 
	head(5) %>% 
	knitr::kable()


# combine with other results results -------------------------------------------

# MSstats intrafraction results:
load(file.path(root,"rdata","msstats_results.rda"))

# annotate with gene symbols
accession <- msstats_results$Protein
idx <- match(accession,gene_map$uniprot)
symbol <- gene_map$symbol[idx]
entrez <- gene_map$entrez[idx]
msstats_results <- tibble::add_column(msstats_results, 
				      Symbol=symbol, .after="Protein")
msstats_results <- tibble::add_column(msstats_results, 
				      Entrez=entrez, .after="Symbol")

# munge harder
contrasts <- msstats_results$Label
msstats_results <- tibble::add_column(msstats_results, 
				      Contrast=contrasts, .before="Protein")
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
final_results = c(results_list[c("F4","F5","F6","F7","F8","F9","F10")],
		  "Control-Mutant"=list(results))
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(final_results,myfile)


## Compare this to our previous results ----------------------------------------

# load the swip data
data(swip_tmt)

# collect edgeR intrafraction results
df0 <- swip_tmt %>% select(Accession,Fraction,PValue,FDR,logFC) %>% unique()
colnames(df0) <- c("protein","contrast","pval","fdr","log2fc")
df0$method="edgeR"

# edgeR summary
message("Summary of edgeR sig prots (FDR<0.05) for intrafraction contrasts:")
df0 %>% group_by(contrast) %>% summarize(n=sum(fdr<0.05)) %>% knitr::kable()

# check number of sig prots:
tmp <- df0 %>% 
	group_by(contrast) %>% 
	filter(fdr<0.05) %>%
       	select(protein,contrast) %>% 
	unique() %>% 
	group_split()
common_prots <- Reduce(intersect,sapply(tmp,"[", "protein"))
all_prots <- Reduce(union,sapply(tmp,"[", "protein"))
message("edgeR: Number of commonly sig DA proteins: ",length(common_prots))
message("edgeR: Total number of unique sig DA proteins: ",length(all_prots))

# edgeR WTvMut
df1 <- swip_tmt %>% 
	select(Accession,Adjusted.PValue,Adjusted.FDR,Adjusted.logFC) %>% 
	unique()
colnames(df1) <- c("protein","pval","fdr","log2fc")
df1$contrast <- "Mutant-Control"
df1$method <- "edgeR"

# edge summary
message("\nSummary of edgeR significant comparisons for Control-Mutant contrast:")
df1 %>% group_by(contrast) %>% 
	summarize(`FDR < 0.05`=sum(fdr<0.05)) %>% knitr::kable()

# MSstats intrafraction
df2 <- msstats_results %>% filter(!is.na(PValue)) %>%
       	select(Protein, Contrast, PValue, log2FC)
colnames(df2) <- c("protein","contrast","pval","log2fc")
df2$contrast <- gsub("Mutant.F[0-9]{1,2}-Control.","",df2$contrast)
df2$method <- "MSstatsTMT"

# edgeR WTvMut
df3 <- results %>% filter(!isSingular) %>% 
	filter(!is.na(PValue)) %>% 
	select(Protein, Contrast, PValue,log2FC)
colnames(df3) <- c("protein","contrast","pval","log2fc")
df3$method <- "MSstatsTMT"

# combine intrafraction restults
df_x = rbind(df0 %>% select(protein,contrast,pval,log2fc,method),df2)

# combine WTvMut results
df_y = rbind(df1 %>% select(protein,contrast,pval,log2fc,method),df3)

# coorelation of pvals
#rho0 <- cor(df_x$pval.edgeR,df_x$pval.MSstats,
#	    method="spearman",use="pairwise.complete")
#rho1 <- cor(df_y$pval.edgeR,df_y$pval.MSstats,
#	    method="spearman",use="pairwise.complete")

#message("Rank-Correlation of PValues:")
#data.table("Control.F#-Mutant.F#" = rho0, "Control-Mutant" = rho1) %>%
#	knitr::kable()

# coorelations of fold change
#rho0 <- cor(df_x$log2fc.edgeR,df_x$log2fc.MSstats,
#	    method="spearman",use="pairwise.complete")
#rho1 <- cor(df_y$log2fc.edgeR,df_y$log2fc.MSstats,
#	    method="spearman",use="pairwise.complete")

#message("Rank-Correlation of log2fc:")
#data.table("Control.F#-Mutant.F#" = rho0, "Control-Mutant" = rho1) %>%
#	knitr::kable()

# the results appear to be highly correlated

## plots ----------------------------------------------------------------------

#df = rbind(df_x %>% select(protein,contrast,pval.edgeR),
#      df_x %>% select(protein,contrast,pval.MSstats))

# How to assess coorelation?
#plot <- ggplot(data=df_x,aes(x=pval.edgeR,y=pval.MSstats))
#plot <- plot + geom_point()
#plot  <- plot + xlab("PValue (edgeR)")
#plot  <- plot + ylab("PValue (MSstatsTMT)")
#plot <- plot + theme(panel.background = element_blank())
#plot <- plot + theme(panel.border = element_rect(colour="black",fill=NA,size=1))

# Plots for intra-fraction contrasts
plot_list <- list()
for (fraction in unique(df_x$contrast)) {
  plot <- ggplot(df_x %>% filter(contrast == fraction), aes(x=pval,colour=method))
  plot <- plot + geom_histogram(bins=100)
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(panel.border = element_rect(colour="black",fill=NA,size=1))
  plot <- plot + ggtitle(fraction)
  plot_list[[fraction]] <- plot
}

# Plot for Control vs Mutant contrast
plot <- ggplot(df_y, aes(x=pval,colour=method))
plot <- plot + geom_histogram(bins=100)
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(colour="black",fill=NA,size=1))
plot <- plot + ggtitle("Control-Mutant")
plot_list[["Control-Mutant"]] <- plot

# save
myfile <- file.path(root,"figs","MSstatsTMT","Pvalue-histograms.pdf")
ggsavePDF(plot_list,myfile)
