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


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- load partition

myfile <- file.path(root,"rdata",input_partition)
df <- data.table::fread(myfile,drop=1)
partition_list <- unlist(apply(df,1,function(x) list(x)),recursive=F,use.names=F)


## ---- save as rda

# add 1 bc python is 0-based 
part <- partition_list[[1]] + 1

# set small modules to 0
modules <- split(names(part),part)
too_small <- as.numeric(names(which(sapply(modules,length) < 5)))
part[part %in% too_small] <- 0

message("Total number of modules: ", length(unique(part))-1)

# save partition
partition <- part
myfile <- file.path(root,"data","partition.rda")
save(partition, file = myfile, version = 2)
