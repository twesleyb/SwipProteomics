#!/usr/bin/env Rscript

# title: Swip Proteomics Plotting
# description: generate module colors
# authors: Tyler W Bradshaw

## ---- INPUT

# input data in root/data
root = "~/projects/SwipProteomics"

part_file = "ne_surprise_partition"

NC_color = "#BEBEBE" # not clustered == "gray"

## ---- OUTPUT 
# * color assignments for every module in graph partition


## ---- FUNCTIONS 

str_to_vec <- function(response) {
  # Parse the python dictionary returned as a string from 
  # system(random_color.py)
	vec <- gsub("'","",gsub("\\[|\\]","",
				trimws(unlist(strsplit(response,",")))))
	return(vec)
}


## ---- Prepare the workspace 

# Load renv
renv::load(root, quiet=TRUE)

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports 
devtools::load_all(root, quiet=TRUE)

# Load TMT data and partition
data(swip)
data(mut_color)
data(list=part_file)


## ---- Generate Module colors 

# The number of colors we need
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
n_colors <- length(modules) # NOTE: M0 will be gray. M(Swip) will be purple.

# Path to python script which is a simple script that uses the python 
# port of randomcolors to generate random colors.
script <- file.path(root,"Py","random_color.py")

# Generate n random colors
cmd <- paste(script,"--count", n_colors, "--luminosity", "bright")
response <- system(cmd, intern = TRUE)

#  Parse the response
colors <- toupper(str_to_vec(response))

if (mut_color %in% colors) { stop("Duplicate colors.") }

# Module color assignments
# Initialize a vector for the module colors
module_colors <- rep(NA,length(modules))
names(module_colors) <- names(modules)

# Insure that M0 is gray and WASH community/module is #B86FAD
module_colors["M0"] <- NC_color
wash_module <- names(which(sapply(modules, function(x) swip %in% x)))
module_colors[wash_module] <- mut_color

# The reamining colors are random
idx <- is.na(module_colors)
module_colors[idx] <- sample(colors,sum(idx))


## ---- Save the data 

message(paste("\nSaving colors."))

# Save updated module colors
namen <- paste0(gsub("partition","colors",part_file),".rda")
myfile <- file.path(root,"data", namen)
save(module_colors,file=myfile,version=2)
