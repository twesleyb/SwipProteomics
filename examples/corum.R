#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab

# prepare renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root,quiet=TRUE)
data(msstats_prot)

# imports
library(geneLists)

data(corum) 

## ---- function 

calcPercentID <- function(pathway) {
  # requires: entrez
  genes <- unlist(pathway,use.names=FALSE)
  percent_id <- sum(genes %in% entrez)/length(genes)
  return(percent_id)
}

# percent identified 
entrez <- unique(msstats_prot$Entrez)
pid <- sapply(corum, calcPercentID)

df <- data.frame("Pathway" = names(corum), 
		 "Size" = sapply(corum,length), 
		 "percentID" = pid)

subdf <- df %>% filter(Size>3,percentID>(2/3))
paths <- unique(subdf$Pathway)

corum_paths <- corum[paths]
entrez <- unlist(corum_paths,use.names=FALSE)
uniprot <- getPPIs::getIDs(entrez,from="entrez",to="uniprot",species="mouse")
corum_prots <- lapply(corum_paths, function(x) as.character(uniprot[x]))

proteins <- unique(msstats_prot$Protein)
mylist <- sapply(corum_prots, function(x) x[x %in% proteins])
idy <- sapply(idx,length)>2
foo = mylist[idy]

data(gene_map)

# fit the module and test the overall contrast
results_list <- list()
for (i in c(1:length(foo))) {
	results_list[[i]] <- fit_module(foo[[i]], msstats_prot, gene_map)
}

results_df$Pathway <- names(foo)

results_df <- do.call(rbind, results_list)
results_df$Padjust <- p.adjust(results_df$Pvalue,method="bonferroni")


fwrite(results_df,"foo.csv")

prots = corum_prots[["CCT complex (chaperonin containing TCP1 complex)"]]
res = fit_module(prots,msstats_prot,gene_map)

data(sig_prots)

sum(prots %in% sig_prots)
