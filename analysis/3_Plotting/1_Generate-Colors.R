#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: generate module colors
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
swip = "Q3UMB9" # uniprot accession of swip.
swip_color = "#B86FAD" # color of swip/wash module.

## OUTPUT:
# * Updated module color assignemnts.

#---------------------------------------------------------------------
## Misc functions
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here)
	}
	root <- here
	return(root)
}

# Parse the python dictionary returned as a string from 
# system(random_color.py)
str_to_vec <- function(response) {
	vec <- gsub("'","",gsub("\\[|\\]","",
				trimws(unlist(strsplit(response,",")))))
	return(vec)
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Load TMT data and partition.
data(tmt_protein)
data(partition)

#---------------------------------------------------------------------
## Generate colors.
#---------------------------------------------------------------------

# The number of colors we need.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
n_colors <- length(modules) - 2 # M0 will be gray. M19 will be purple.

# Path to python script which is a simple script that uses the python 
# port of randomcolors to generate random colors.
script <- file.path(root,"Py","random_color.py")

# Generate n random colors.
cmd <- paste(script,"--count", n_colors)
response <- system(cmd, intern = TRUE)

#  Parse the response.
colors <- toupper(str_to_vec(response))

if (swip_color %in% colors) { stop("Duplicate colors.") }

# Module color assignments.
# Initialize a vector for the module colors.
module_colors <- rep(NA,length(modules))
names(module_colors) <- names(modules)

# Insure that M0 is gray and WASH community/module is #B86FAD.
module_colors["M0"] <- col2hex("gray")
wash_module <- names(which(sapply(modules, function(x) swip %in% x)))
module_colors[wash_module] <- swip_color

# The reamining colors are random.
idx <- is.na(module_colors)
module_colors[idx] <- sample(colors,sum(idx))

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

message(paste("\nSaving colors."))

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
