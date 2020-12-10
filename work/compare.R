#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(swip_partition)
data(swip_sig_prots); s <- sig_prots
data(msstats_sig_prots); m <- sig_prots

library(dplyr)
library(data.table)

S <- split(names(partition) %in%  s, partition)
M <- split(names(partition) %in%  m, partition)
names(M) <- names(S) <- paste0("M",names(M))

sapply(M, sum) %>% knitr::kable()

sapply(S, sum) %>% knitr::kable()

modules <- split(names(partition),partition)[-1]

data(msstats_gene_map)

results <- list()
for (i in c(1:length(modules))){
  prots <- modules[[i]]
  u <- intersect(prots[prots %in% s],prots[prots %in% m])
  x = prots[prots %in% s]
  to_rm <- x[x %notin% m] # to remove
  if (length(to_rm) > 0) {
	  to_rm = mapID(to_rm,"uniprot","symbol")
  }
  y = prots[prots %in% m]
  to_add <- y[y %notin% s] # to add
  if (length(to_add) > 0)  {
	  to_add = mapID(to_add,"uniprot","symbol")
  }
  df <- data.table(current=sum(prots %in% s), to_label = sum(prots %in% m), to_rm = length(to_rm), to_add = length(to_add))
  results[[i]] <- list(summary = df, to_rm = to_rm, to_add = to_add)
}

results
