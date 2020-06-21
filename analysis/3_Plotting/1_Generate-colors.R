#!/usr/bin/env Rscript

## OPTIONS:
min_size = 5 # minimum size of a module.
max_size = 100 # maximum size of a module.
seed = 7 # Seed for reproducibility.
swip = "Q3UMB9" # uniprot accession of swip.

## OUTPUT:
# * Updated module color assignemnts.

# To conserve on colors. All communities will be the same color.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

start <- Sys.time()
message(paste("Starting analysis at:", start))

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(igraph)
	library(WGCNA)
	library(colorspace)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Load TMT data and partition.
data(tmt_protein)
data(partition)

set.seed(seed)

#---------------------------------------------------------------------
## Load the initial partition of the graph. 
#---------------------------------------------------------------------
# Clean-up intial partition of the graph into large communties
# (>100 nodes).

# All modules.
modules <- split(names(partition),partition)

# All communities. 
# Load the intial partition of the graph into large communities.
myfile <- file.path(root,"rdata","Swip_initial_partition.csv")
cpartition <- fread(myfile,drop=1) %>% unlist() + 1 # Add 1 because python uses 0 index.
communities <- split(names(cpartition),cpartition)
names(communities) <- paste0("C",names(communities))

# cmodules is the modules contained by every community.
cmodules <- lapply(communities, function(x) paste0("M",unique(partition[x])))

# combine modules that only contain "M0"
idx <- as.numeric(which(sapply(cmodules, function(x) all(x == "M0"))))
new_part <- cpartition
new_part[new_part %in% idx] <- 0
new_cpart <- reset_index(new_part)
cpartition <- new_cpart

# New communtities
communities <- split(names(new_cpart),new_cpart)
names(communities) <- paste0("C",names(communities))

# Save as rda.
save(communities,file=file.path(root,"data","communities.rda"),version=2)
save(cpartition,file=file.path(root,"data","cpartition.rda"),version=2)

#---------------------------------------------------------------------
## Cluster the communities.
#---------------------------------------------------------------------

#---------------------------------------------------------------------
## Generate colors.
#---------------------------------------------------------------------
# Modules will be assigned the same color as their community.

message(paste("\nColors were generated online at 'coolars.co'"))
data(coolors)

# Generate community colors.
n_comm <- length(communities) -1
community_colors <- c(col2hex("gray"), sample(coolors,n_comm))
names(community_colors) <- names(communities)

# Insure that WASH community/module is #B86FAD
wash_community <- names(which(sapply(communities, function(x) swip %in% x)))
community_colors[wash_community] <- "#B86FAD"

# Protein color assignments.
protein_colors <- community_colors[paste0("C",cpartition)]
names(protein_colors) <- names(cpartition)

# Module color assignments.
module_colors <- sapply(split(protein_colors,
			      partition[names(protein_colors)]),unique)
names(module_colors) <- paste0("M",names(module_colors))
module_colors <- unlist(module_colors,use.names=TRUE,recursive=FALSE)
module_colors <- unlist(module_colors,use.names=FALSE) 
names(module_colors) <- names(module_colors)

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# Save updated community colors.
myfile <- file.path(root,"data","community_colors.rda")
save(community_colors,file=myfile,version=2)

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
