#!/usr/bin/env Rscript

# title: SwipProteomics
# description:
# author: twab

## ---- INPUTs 
input_part = "ne_surprise_partition"
input_colors = "ne_surprise_colors"

# * msstats_prot
# * partition
# * colors


## ---- OPTIONs 
fig_width = 5
fig_height = 5


## ---- OUTPUTs 
# * pdf of protein pca plot with module colors


## ---- Prepare the workspace 

# Load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Local Imports
devtools::load_all(root, quiet=TRUE)

# Project directories:
figsdir <- file.path(root,"figs","Proteins")

if (!dir.exists(figsdir)) { 
  dir.create(figsdir) 
  message("mkdir ", figsdir) 
}


## ---- Prepare the data for ploting 

# load partition
data(list=input_part)

# load the data
data(msstats_prot)

# load colors
data(list=input_colors)

# all modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# coerce tidy data to a matrix
# summarize three replicates as median
dm <- msstats_prot %>% 
	filter(Protein %in% names(partition)) %>%
	group_by(Protein, Condition) %>% 
	summarize(med_Abundance = median(Abundance),.groups="drop") %>%
	reshape2::dcast(Protein ~ Condition, value.var= "med_Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

# normalize each protein such that its sum is 1.
norm_dm <- t(apply(dm, 1, function(x) x/sum(x)))

# Drop un-clustered proteins
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
# Collect data for plotting
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")
# Generate the plot
plot <- ggplot(df)
plot <- plot + aes(x, y)
plot <- plot + geom_point()
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
# Get plot's data, annotate with protein module membership
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
