#!/usr/bin/env Rscript 

# title: SwipProteomics
# description: generate protein co-variation (correlation) network
# author: twab

## ---- Input:
root <- "~/projects/SwipProteomics"

## ---- Output:

# * adjm.rda
# * ne_adjm.rda

# * ne_adjm.csv --> for leidenalg clustering!

# NOTE: large ouput files saved in root/rdata bc too big to be tracked by git


## ---- prepare the working environment

# load renv
renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load data in root/data
data(gene_map)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(neten) # twesleyb/neten
  library(igraph)
  library(data.table)
})


## ---- create covariation network

message("Generating covariation network...")

# summarize median of three mixtures
dm <- msstats_prot %>%  
	group_by(Protein, Condition) %>%
	summarize(med_Abundance = median(Abundance),
		  .groups="drop") %>% 
	reshape2::dcast(Protein ~ Condition, value.var = "med_Abundance") %>% 
	as.data.table() %>%
	as.matrix(rownames="Protein")


# there are a small number proteins with some missing vals 
# e.g. Q9QUN7 = low abundance only quantified in 4/7 fractions
# remove these proteins
idx <- apply(dm,1,function(x) any(is.na(x)))
warning(sum(idx)," proteins with any missing values are removed.")
filt_dm <- dm[!idx,]


## calculate coorrelation matrix
adjm <- cor(t(filt_dm), method="pearson",use="complete.obs")


## ---- network enhancement
# REF: Wang et al., 2018 (Nature Communications)

message("\nPerforming network enhancement...")

ne_adjm <- neten(adjm) # result is robust to neten parameters


## ---- save network

# coerce to data.table and write to file
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","adjm.csv")
data.table::fwrite(adjm_dt, myfile)

# coerce to data.table and write to file
ne_adjm_dt <- as.data.table(ne_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ne_adjm.csv")
data.table::fwrite(ne_adjm_dt, myfile)


## ---- save as rda

# adjm
myfile <- file.path(root,"rdata","adjm.rda")
save(adjm, file=myfile,version=2)

# ne adjm
myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm, file=myfile,version=2)
