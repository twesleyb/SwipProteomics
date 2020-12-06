#!/usr/bin/env bash

# load renv, project functions and data
renv()

data(gene_map)
data(lopit_prots)
data(ne_surprise_partition)


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






