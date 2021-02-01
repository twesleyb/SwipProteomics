#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: generate WT and MUT protein co-variation (correlation) network 
# and perform network enhancement in preparation for clustering and perm testing

## ---- Input:
root <- "~/projects/SwipProteomics"
devtools::load_all(root)

## ---- Output:

# * wt_adjm.rda
# * wt_ne_adjm.rda
# * mut_adjm.rda
# * mut_ne_adjm.rda

# * wt_ne_adjm.csv --> for leidenalg clustering!
# * mut_ne_adjm.csv --> for leidenalg clustering!

# NOTE: large files are saved in root/rdata bc too big to be tracked by git


## ---- prepare the working environment

# load data in root/data
data(swip_tmt)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(neten) # twesleyb/neten
  library(igraph)
  library(data.table)
})


## ---- create covariation network

message("Generating covariation network...")

# collect control (WT) data as matrix
wt_dm <- swip_tmt %>% filter(grepl("Control",Condition)) %>%
	reshape2::dcast(Protein ~ Mixture + Condition, 
			value.var = "Abundance") %>%
	as.data.table() %>%
	as.matrix(rownames="Protein")

# collect WASHC4 mutant data as matrix
mut_dm <- swip_tmt %>% filter(grepl("Mutant",Condition)) %>%
	reshape2::dcast(Protein ~ Mixture + Condition, 
			value.var = "Abundance") %>%
	as.data.table() %>%
	as.matrix(rownames="Protein")

# there should be no missing values
stopifnot(!any(apply(wt_dm, 1, function(x) any(is.na(x)))))
stopifnot(!any(apply(mut_dm, 1, function(x) any(is.na(x)))))

# calculate correlation matrix
wt_adjm <- cor(t(wt_dm), method="pearson",use="complete.obs")
mut_adjm <- cor(t(mut_dm), method="pearson",use="complete.obs")


## ---- perform network enhancement

# Wang et al., 2018 (Nature Communications; PMID:30082777)
# NOTE: result is robust to neten parameters

message("Performing network enhancement...")
wt_netw <- neten(wt_adjm) 

message("Performing network enhancement...")
mut_netw <- neten(mut_adjm)


## ---- save networks as csv in rdata

# save WT adjm
myfile <- file.path(root,"rdata","wt_adjm.csv")
wt_adjm %>% as.data.table(keep.rownames="Protein") %>% 
	data.table::fwrite(myfile)
message("saved: ", myfile)

# save MUT adjm
myfile <- file.path(root,"rdata","mut_adjm.csv")
mut_adjm %>% as.data.table(keep.rownames="Protein") %>% 
	data.table::fwrite(myfile)
message("saved: ", myfile)

# save WT netw
myfile <- file.path(root,"rdata","wt_netw.csv")
wt_netw %>% as.data.table(keep.rownames="Protein") %>%
	data.table::fwrite(myfile)
message("saved: ", myfile)

# save MUT netw
myfile <- file.path(root,"rdata","mut_netw.csv")
mut_netw %>% as.data.table(keep.rownames="Protein") %>%
	data.table::fwrite(myfile)
message("saved: ", myfile)


## ---- save data as rda

# save wt adjm
myfile <- file.path(root,"rdata","wt_adjm.rda")
save(wt_adjm, file=myfile,version=2)
message("saved: ", myfile)

# save mut adjm
myfile <- file.path(root,"rdata","mut_adjm.rda")
save(mut_adjm, file=myfile,version=2)
message("saved: ", myfile)

# save WT netw
myfile <- file.path(root,"rdata","wt_netw.rda")
save(wt_netw, file=myfile,version=2)
message("saved: ", myfile)

# save MUT netw
myfile <- file.path(root,"rdata","mut_netw.rda")
save(mut_netw, file=myfile,version=2)
message("saved: ", myfile)
