#!/usr/bin/env Rscript

## INPUTs ----------------------------------------------------------------------
# * msstats_prot
# * partition
# * colors

# OPTIONs ---------------------------------------------------------------------
fig_width = 5
fig_height = 5

## OUTPUTs --------------------------------------------------------------------
# * pdf of protein pca plot with module colors.

# FUNCTIONs -------------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# Get the repository's root directory.
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

## Prepare the workspace ------------------------------------------------------

# Load renv
root <- getrd()
renv::load(root,quiet=TRUE)

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Local Imports
suppressMessages({ devtools::load_all() })

# Project directories:
figsdir <- file.path(root,"figs","Proteins")


## Prepare the data for ploting ------------------------------------------------

# load partition
data(partition)

# load the data
data(msstats_prot)

# load colors
data(module_colors)

# all modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# Drop QC and coerce data to data.matrix.
dm <- msstats_prot %>% 
	filter(Protein %in% names(partition)) %>%
	as.data.table() %>%
	reshape2::dcast(Protein ~ interaction(Mixture, BioFraction, 
					      Condition), value.var= "Abundance") %>%
		as.data.table() %>% as.matrix(rownames="Protein")

# Scale rows (proteins).
# NOTE: normalization is applied row-wise (dim=1), but we
# need to transpose the output such that it is the same
# dimensions as the input.
norm_dm <- t(apply(dm,1,function(x) x/sum(x)))
#norm_dm = dm

# Drop un-clustered proteins.
idx <- rownames(norm_dm) %in% modules[["M0"]]
filt_dm <- norm_dm[!idx,]

# do pca
pca <- prcomp(filt_dm)
pca_summary <- as.data.frame(t(summary(pca)$importance))

# get top 2
idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
pca_summary <- pca_summary[idx,]
top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
names(top2_pc) <- head(rownames(pca_summary),2)

# Plot axis labels:
x_label <- paste0(names(top2_pc)[1],
		  " (PVE: ",round(100*top2_pc[1],2)," %)")
y_label <- paste0(names(top2_pc)[2],
		  " (PVE: ",round(100*top2_pc[2],2)," %)")

# Collect data for plotting.
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")

# Generate the plot.
plot <- ggplot(df, aes(x,y)) + geom_point()
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))

# Get plot's data, annotate with protein module membership.
df <- plot$data
df$module <- paste0("M",partition[rownames(df)])

# Annotate with color assignment
df$color <- module_colors[df$module]

# Add color to plot
plot <- plot + geom_point(data = df, aes(colour=factor(module)))
plot <- plot + scale_colour_manual(values=df$color)

# Drop the legend
plot <- plot + theme(legend.position = "none")

# Save to file
myfile <- file.path(figsdir,"Protein_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
