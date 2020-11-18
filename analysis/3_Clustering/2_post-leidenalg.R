#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## input
input_partition <- "ne_cpm_surprise_partition.csv"

## output
output_partition <- "leidenalg_partition.rda"


## ---- prepare the env

root <- "~/projects/SwipProteomics"
renv::load(root)

# library(SwipProteomics)
#devtools::load_all(root)


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- load partition

myfile <- file.path(root,"rdata",input_partition)
partition <- data.table::fread(myfile,drop=1) %>% unlist()


## ---- save

myfile <- file.path(root,"data",output_partition)
save(partition, file = myfile, version = 2)
