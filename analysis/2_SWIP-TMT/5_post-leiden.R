#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## inputs:
root <- "~/projects/SwipProteomics"
part_file <- file.path(root, "rdata", "leidenalg_partition.csv")

## output:
output_partition <- file.path(root, "data", "leidenalg_partition.rda")

## prepare env
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# load partition and save as named vector
part <- data.table::fread(part_file, drop = 1)
partition <- unlist(part) # coerce to named numeric
save(partition, file = output_partition, version = 2)
