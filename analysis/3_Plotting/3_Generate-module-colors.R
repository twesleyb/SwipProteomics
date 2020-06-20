#!/usr/bin/env Rscript

## OPTIONS:
min_size = 5

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
	library(colorspace)
	library(data.table)
	library(dplyr)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Load TMT data and partition.
data(tmt_protein)
data(partition)

#---------------------------------------------------------------------
## Collect modules and communties.
#---------------------------------------------------------------------

# All communities. 
# Load the intial partition of the graph into large communities.
myfile <- file.path(root,"rdata","Swip_initial_partition.csv")
community_part <- fread(myfile,drop=1) %>% unlist()
communities <- split(names(community_part),community_part)
names(communities) <- paste0("C",names(communities))

# Remove communties that are smaller than min_size.
idx <- sapply(communities,length) < min_size
names(communities)[idx] <- "C0"

# Collect unclustered proteins, and reset their Community membership to 0.
not_clustered <- unlist(communities[names(communities) == "C0"])
community_part[not_clustered] <- 0

# Reset partition index.
new_part <- reset_index(community_part)
communities <- split(names(new_part),new_part)
names(communities) <- paste0("C",names(communities))

# Sizes of the communities.
community_sizes <- sapply(communities,length)
to_split <- names(which(community_sizes > 100))
to_split <- to_split[!to_split=="C0"]
n_split <- length(to_split)

# Proteins assigned to each module.
community_prots <- split(partition,community_part)

# All modules.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

#---------------------------------------------------------------------
## Generate colors.
#---------------------------------------------------------------------

# Generate community colors.
n_communities <- length(communities) -1
community_colors <- c(col2hex("gray"),colorspace::rainbow_hcl(n_communities))
names(community_colors) <- names(communities)

# Generate module colors.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))
n_modules <- length(modules) -1
all_colors <- c(col2hex("gray"),colorspace::rainbow_hcl(n_modules))

# Insure that WASH module is #B86FAD
#swip = "Q3UMB9"
#m <- paste0("M",partition[swip])
#module_colors[m] <- "#B86FAD"

#--------------------------------------------------------------------
## Organize the colors.
#--------------------------------------------------------------------

#  Calculate module eigengenes.
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()
ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Extract ME matrix.
ME_dm <- as.matrix(ME_data$eigengenes)

# Compute relationships between modules.
ME_adjm <- WGCNA::bicor(ME_dm)

# Cluster heirarchically.
d <- as.dist(1 - ME_adjm)
hc <- hclust(d)

# Assign names to colors in heirarchical order.
names(all_colors) <- hc$labels[hc$order]

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# Save updated community colors.
myfile <- file.path(root,"data","community_colors.rda")
save(community_colors,file=myfile,version=2)

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)
