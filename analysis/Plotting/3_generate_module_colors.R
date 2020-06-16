#!/usr/bin/env Rscript

# Colors were generated online at: https://coolors.co/.

## OUTPUT:
# * Updated module color assignemnts.

## R Options:
options(renv.config.synchronized.check = FALSE) # skip renv::check(repo).
options(renv.settings.snapshot.type = "simple") # use simple renv::snapshot.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(colorspace)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Load partition.
data(partition)

# Generate colors.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))
n_modules <- length(modules) -1
module_colors <- c(col2hex("gray"),colorspace::rainbow_hcl(n_modules))
names(module_colors) <- names(modules)


# Insure that WASH module is #B86FAD
swip = "Q3UMB9"
m <- paste0("M",partition[swip])
module_colors[m] <- "#B86FAD"

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
