#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## input in root/rdata
input_partition <- "ne_surprise_surprise_partition.csv"
# output is [input_partition].rda


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
df <- data.table::fread(myfile,drop=1)
partition_list <- unlist(apply(df,1,function(x) list(x)),recursive=F,use.names=F)


## ---- save

partition <- partition_list[[1]]

myfile <- paste0(tools::file_path_sans_ext(input_partition),".rda")
output_file <- file.path(root,"data",myfile)
save(partition, file = output_file, version = 2)
