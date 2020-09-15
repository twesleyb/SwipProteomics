#!/usr/bin/env Rscript

# OPTIONS:
fig_width = 5
fig_height = 5

## OUTPUT:
# * pdf of module protein pca plot.

## Misc functions -------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# get a git repo's root directory
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

# global imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# local imports
suppressMessages({ devtools::load_all() })

# project directories:
figsdir <- file.path(root,"figs","Samples")

# if necessary, create figsdir
if (!dir.exists(figsdir)) {
	dir.create(figsdir)
}

## Prepare the data for ploting -----------------------------------------------

# load the proteomics data
data(swip_tmt)

# cast into a matrix
dm <- tmt_protein %>%  
	reshape2::dcast(Accession ~ Sample, value.var= "Intensity") %>%
		as.data.table() %>% as.matrix(rownames="Accession")

# PCA
pca <- prcomp(log2(t(dm))) # log2
pca_summary <- as.data.frame(t(summary(pca)$importance))
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

# annotate with group info
idx <- match(rownames(df),tmt_protein$Sample)
df$group <- interaction(tmt_protein$Genotype[idx],tmt_protein$Fraction[idx])

# generate the plot
plot <- ggplot(df, aes(x,y,color=group)) + geom_point()
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

# Save to file.
myfile <- file.path(figsdir,"Sample_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
