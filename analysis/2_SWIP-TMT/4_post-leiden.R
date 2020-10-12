#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## inputs:

root <- "~/projects/SwipProteomics"

# input partitions in root/rdata
part_file1 <- file.path(root,"rdata","leidenalg_partition.csv")
part_file2 <- file.path(root,"rdata","multiplex_partition.csv")

## output:

# partition R objects saved in root/data
part_out1 <- file.path(root,"data","leidenalg_partition.rda")
part_out2 <- file.path(root,"data","multiplex_partition.rda")


## prepare env
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# local imports
#devtools::load_all(quiet=TRUE)

# the msstats protein-level data
#data(msstats_prot)

# load partition and save as named vector
part <- data.table::fread(part_file1,drop=1)
partition <- unlist(part) # coerce to named numeric
save(partition,file=part_out1,version=2)

# repeat for multiplex partition
part <- data.table::fread(part_file2,drop=1)
partition <- unlist(part) # coerce to named numeric
save(partition,file=part_out2,version=2)
