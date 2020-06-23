#!/usr/bin/env Rscript

## OPTIONS:
seed = 7 # Seed for reproducibility.
swip = "Q3UMB9" # uniprot accession of swip.
swip_color = "#B86FAD" # color of wash module.

## OUTPUT:
# * Updated module color assignemnts.

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

data(coolors)
data(colormind)

#---------------------------------------------------------------------
## Generate colors.
#---------------------------------------------------------------------

# All modules.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# Module color assignments -- combine coolors and colormind.
colors <- c(coolors, sample(colormind,length(modules)-length(coolors)))
module_colors <- sample(colors,length(modules))
names(module_colors) <- names(modules)

# Insure that WASH community/module is #B86FAD
wash_module <- names(which(sapply(modules, function(x) swip %in% x)))
module_colors[wash_module] <- swip_color

# Insure that M0 is gray.
module_colors["M0"] <- col2hex("gray")

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
