#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up partitions from running leidenalg

## input:

root <- "~/projects/SwipProteomics"

part_file1 <- file.path(root,"rdata","leidenalg_partition.csv")
part_file2 <- file.path(root,"rdata","multi_leidenalg_partition.csv")

## output:

part_out1 <- file.path(root,"rdata","leidenalg_partition.rda")
part_out2 <- file.path(root,"rdata","multi_leidenalg_partition.rda")

# * norm_prot.rda in root/data

## prepare env
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# local imports
devtools::load_all(quiet=TRUE)

# swip's uniprot
data(swip)

# the msstats protein-level data
data(msstats_prot)

# load partition and save as named vector
part <- data.table::fread(part_file1,drop=1)
partition <- unlist(part) # coerce to named numeric
save(partition,file=part_out1,version=2)

# repeat for multiplex partition
part <- data.table::fread(part_file2,drop=1)
partition <- unlist(part) # coerce to named numeric
save(partition,file=part_out2,version=2)


## create norm_prot for permutation testing -----------------------------------

# load the data and save as norm protein for permutation testing
norm_prot <- msstats_prot %>% as.data.table() %>%
	dcast(Protein ~ interaction(Mixture,Channel,Genotype),value.var="Abundance") %>%
	na.omit() %>% as.matrix(rownames="Protein")

# save
myfile <- file.path(root,"rdata","msstats_norm_prot.rda")
save(norm_prot,file=myfile,version=2)
