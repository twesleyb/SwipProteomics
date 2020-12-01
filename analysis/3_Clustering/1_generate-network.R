#!/usr/bin/env Rscript 

# title: SwipProteomics
# description: generate protein co-variation (correlation) network
# author: twab

## ---- Input:
root <- "~/projects/SwipProteomics"

# * msstats_prot and other R data in root/data

## Options:
# rm proteins with poor fit?
rm_poor <- TRUE 

# network enhancement params
alpha_param = 0.9
diffusion_param = 1.0


## ---- Output:
# adjm.rda
# ppi_adjm.rda
# ne_adjm.rda
# ne_adjm.csv --> for leidenalg clustering!

# NOTE: large ouput files saved in root/rdata bc too big to be tracked by git


## ---- prepare the working environment

# load renv
renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load data in root/data
data(gene_map)
data(poor_prots)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(neten) # twesleyb/neten
  library(igraph)
  library(getPPIs) # twesleyb/getPPIs
  library(data.table)
})


# load mouse PPIs compiled from HitPredict
data(musInteractome) 


## ---- create covariation network

dm <- msstats_prot %>% 
	mutate(Intensity = 2^Abundance) %>% 
	group_by(Protein) %>%
	mutate(scale_Intensity = Intensity/sum(Intensity)) %>%
	group_by(Protein, Condition) %>%
	summarize(med_Intensity = log2(median(scale_Intensity)), .groups="drop") %>% 
	reshape2::dcast(Protein ~ Condition, value.var = "med_Intensity") %>% 
	as.data.table() %>%
	as.matrix(rownames="Protein")

# there are a small number proteins with some missing vals 
# e.g. Q9QUN7 = low abundance only quantified in 4/7 fractions
# remove these proteins
idx1 <- apply(dm,1,function(x) any(is.na(x)))
warning(sum(idx1)," proteins with any missing values are removed.")
filt_dm <- dm[!idx1,]

## calculate coorrelation matrix
adjm <- cor(t(filt_dm), method="pearson",use="complete.obs")


## ---- network enhancement
# REF: Wang et al., 2018 (Nature Communications)

ne_adjm <- neten(adjm, alpha = alpha_param, diffusion = diffusion_param)


## ---- create ppi network

# NOTE: PPIs are NOT used to identify communities

# map uniprot to entrez
uniprot <- colnames(adjm)
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]

# given entrez, collect ppis from musInteractome
ppi_df <- musInteractome %>% 
	subset(osEntrezA %in% entrez & osEntrezB %in% entrez) %>% 
	select(osEntrezA, osEntrezB)

# map back to uniprot and cast to matrix for igraph
idx <- match(ppi_df$osEntrezA,gene_map$entrez)
idy <- match(ppi_df$osEntrezB,gene_map$entrez)
ppi_dm <- ppi_df %>% 
	mutate(ProtA = gene_map$uniprot[idx], ProtB = gene_map$uniprot[idy]) %>%
	select(ProtA,ProtB) %>% as.matrix()

# create igraph graph
g <- igraph::graph_from_edgelist(ppi_dm, directed=FALSE)

# simplify (weight=0,1) and get the adjacency matrix
ppi_adjm <- as.matrix(igraph::as_adjacency_matrix(igraph::simplify(g)))

# collect proteins that are missing
missing_prots <- colnames(adjm)[colnames(adjm) %notin% colnames(ppi_adjm)]
# these proteins are unconnected^, but we include them in the ppi_adjm so the
# networks have matching vertex sets

# add missing cols
tmp_cols <- matrix(0, nrow=nrow(ppi_adjm),ncol=length(missing_prots))
colnames(tmp_cols) <- missing_prots
rownames(tmp_cols) <- rownames(ppi_adjm)
tmp_dm <- cbind(ppi_adjm,tmp_cols)

# add missing rows
tmp_rows <- matrix(0, nrow=length(missing_prots),ncol=ncol(tmp_dm))
colnames(tmp_rows) <- colnames(tmp_dm)
rownames(tmp_rows) <- missing_prots

# full(ppi)_adjm
full_adjm <- rbind(tmp_dm,tmp_rows)

# sort rows and cols to match adjm
ppi_adjm <- full_adjm[colnames(adjm),rownames(adjm)]


## ---- save networks

# coerce to data.table and write to file
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","adjm.csv")
data.table::fwrite(adjm_dt, myfile)

# coerce to data.table and write to file
ne_adjm_dt <- as.data.table(ne_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ne_adjm.csv")
data.table::fwrite(ne_adjm_dt, myfile)

# coerce to data.table and write to file
ppi_dt <- as.data.table(ppi_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ppi_adjm.csv")
data.table::fwrite(ppi_dt, myfile)


## ---- save as rda

# adjm
myfile <- file.path(root,"rdata","adjm.rda")
save(adjm, file=myfile,version=2)

# ne adjm
myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm, file=myfile,version=2)

# ppi adjm
myfile <- file.path(root,"rdata","ppi_adjm.rda")
save(ppi_adjm, file=myfile,version=2)
