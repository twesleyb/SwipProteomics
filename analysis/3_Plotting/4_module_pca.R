#!/usr/bin/env Rscript

# Colors were generated online at: https://coolors.co/.

# OPTIONS:
fig_width = 5
fig_height = 5

## OUTPUT:
# * pdf of module protein pca plot.

## R Options:
options(renv.config.synchronized.check = FALSE) # skip renv::check(repo).
options(renv.settings.snapshot.type = "simple") # use simple renv::snapshot.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

# Project directories:
figsdir <- file.path(root,"figs","Modules")

# If necessary, create figsdir.
if (!dir.exists(figsdir)) {
	dir.create(figsdir)
}

# Load partition.
data(partition)

# Load colors.
data(module_colors)

#---------------------------------------------------------------------
## Prepare the data for ploting.
#---------------------------------------------------------------------

# All modules.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# Load the proteomics data.
data(tmt_protein)

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

# Get plot's data, annotate with protein module membership.
df <- plot$data
df$module <- paste0("M",partition[rownames(df)])

# Annotate with color assignment.
df$color <- module_colors[df$module]

# Add color to plot.
plot <- plot + geom_point(data = df, aes(colour=factor(module)))
plot <- plot + scale_colour_manual(values=df$color)

# Drop the legend.
plot <- plot + theme(legend.position = "none")

# Save to file.
myfile <- file.path(figsdir,"Module_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
