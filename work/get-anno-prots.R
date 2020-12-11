#!/usr/bin/env bash

# load renv, project functions and data
renv()

data(gene_lists)
data(swip_gene_map)
data(swip_partition)

namen <- names(gene_lists)

path_entrez = gene_lists[["LopitDC: PM"]]
path_uniprot = getPPIs::getIDs(path_entrez,"entrez","uniprot","mouse")
length(path_uniprot)
length(path_entrez)

module = names(which(partition==4))
sum(module %in% path_uniprot)

x = gene_lists


myfun <- function(path,m){
  pathway = lopit_prots[[toupper(path)]]
  if (is.null(pathway)) { 
	stop("path must be one of:\n",
	     paste(names(lopit_prots),collapse=", ")) 
  }
  prots = partition[partition==m] %>% names()
  idx = prots %in% pathway
  stopifnot(any(idx))
  mapID(prots[idx],"uniprot","symbol")
}

myfun("lysosome", 21)


data(corum_prots)



data(lysosome)

lysosome

df <- data.table(Protein = lysosome, anno = "LopitDC:Lysosome")
fwrite(df,"foo.csv")






