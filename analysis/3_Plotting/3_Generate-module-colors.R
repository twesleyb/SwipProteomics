#!/usr/bin/env Rscript

## OPTIONS:
min_size = 5

## OUTPUT:
# * Updated module color assignemnts.

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
module_colors <- c(col2hex("gray"),colorspace::rainbow_hcl(n_modules))
names(module_colors) <- names(modules)

# Insure that WASH module is #B86FAD
swip = "Q3UMB9"
m <- paste0("M",partition[swip])
module_colors[m] <- "#B86FAD"

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# Save updated community colors.
myfile <- file.path(root,"data","community_colors.rda")
save(community_colors,file=myfile,version=2)

# Save updated module colors.
myfile <- file.path(root,"data","module_colors.rda")
save(module_colors,file=myfile,version=2)

#--------------------------------------------------------------------
## Generate a plot.
#--------------------------------------------------------------------

# Working with tmt_protein...
# Drop QC and coerce data to data matrix.
dm <- tmt_protein %>% filter(Treatment != "SPQC") %>% 
	filter(Accession %in% names(partition)) %>%
	as.data.table() %>%
	reshape2::dcast(Accession ~ Sample, value.var= "Intensity") %>%
		as.data.table() %>% as.matrix(rownames="Accession")

# Scale rows.
# NOTE: normalization is applied row-wise (dim=1), but we
# need to transpose the output such that it is the same
# dimensions as the input.
norm_dm <- t(apply(dm,1,function(x) x/sum(x)))

# Drop un-clustered proteins.
idx <- rownames(norm_dm) %in% modules[["M0"]]
norm_dm <- norm_dm[!idx,]

# Generate plot.
plot <- ggplotPCAprot(norm_dm,scale=TRUE,center=TRUE)

# Collect plots data and annotate with color assignments.
df <- plot$data
df$community <- paste0("C",community_part[rownames(df)])
df$module <- paste0("M",partition[rownames(df)])
df$community_color <- community_colors[df$community]
df$module_color <- module_colors[df$module]

#--------------------------------------------------------------------
## How to organize the colors?
#--------------------------------------------------------------------

# 57 communities.
# 250 modules.

plot <- plot + geom_point(data=df, aes(colour=factor(community)))
plot <- plot + scale_colour_manual(values=df$community_color)
plot <- plot + theme(legend.position = "none")
plot

