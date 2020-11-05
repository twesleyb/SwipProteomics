#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable
# and saves the results in root/data as R object

## inputs:
root <- "~/projects/SwipProteomics"
#part_file <- file.path(root, "rdata", "leidenalg_partition.csv")
part_file <- file.path(root,"rdata", "ne_rber_surprise_partition.csv")
resolution <- 1

## output:
output_partition <- file.path(root, "data", "leidenalg_partition.rda")

## prepare env
renv::load(root)
devtools::load_all(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# load partition and save as named vector
part <- data.table::fread(part_file, drop = 1)
partition <- unlist(part[resolution,]) # coerce to named numeric
save(partition, file = output_partition, version = 2)
