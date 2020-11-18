#!/usr/bin/env Rscript

# OPTIONS:
fig_width = 5
fig_height = 5

## OUTPUT ---------------------------------------------------------------------
# * pdf of module protein pca plot in figs/Modules

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
data(msstats_prot)

# cast into a matrix
dm <- msstats_prot %>%
	reshape2::dcast(Protein ~ interaction(Mixture, BioFraction, Genotype), 
			value.var= "norm_Abundance") %>%
		as.data.table() %>% as.matrix(rownames="Protein")

# there should be no missing values
to_drop <- apply(dm,1,function(x) any(is.na(x)))
subdm <- dm[!to_drop,]

# PCA
pca <- prcomp(t(subdm)) # log2
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
geno <- sapply(strsplit(rownames(df),"\\."),"[",3)
frac <- sapply(strsplit(rownames(df),"\\."),"[",2)
df$group <- interaction(geno,frac)

# generate the plot
plot <- ggplot(df, aes(x,y,color=group)) + geom_point(size=4)
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

###############################################################################
## repeate with raw data -- show batch effect

myfile <- file.path(root,"rdata","msstats_psm.rda")
load(myfile)

# cast into a matrix
dm <- msstats_psm %>%
	reshape2::dcast(PSM ~ Mixture + Channel + Condition, 
			value.var= "Intensity") %>%
		as.data.table() %>% as.matrix(rownames="PSM")

# there should be no missing values
to_drop <- apply(dm,1,function(x) any(is.na(x)))
subdm <- dm[!to_drop,]

# PCA
pca <- prcomp(t(log2(subdm))) # log2
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
condition <- sapply(strsplit(rownames(df),"\\_"),"[",3)
df$group <- condition

# generate the plot
plot <- ggplot(df, aes(x,y,color=group)) + geom_point(size=4)
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


# Save to file.
myfile <- file.path(figsdir,"raw_Sample_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)

## ------------------------------------------------------------------

## SPN normalization for MSstatsProt
#data(swip)
#df <- msstats_prot %>% filter(Protein==swip) 


quit()

## ------------------------------------------------------------------
## what does old norm data look like?
## ~90% PVE >> msstats --> need to include mixture/batch normalization?
data(swip_tmt)

# cast into a matrix
dm <- swip_tmt %>%
	reshape2::dcast(Accession ~ interaction(Experiment, Fraction, Genotype), 
			value.var= "Intensity") %>%
		as.data.table() %>% as.matrix(rownames="Accession")
to_drop <- apply(dm,1,function(x) any(is.na(x)))
subdm <- dm[!to_drop,]

# PCA
pca <- prcomp(t(log2(subdm))) # log2 if necessary
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
geno <- sapply(strsplit(rownames(df),"\\."),"[",3)
frac <- sapply(strsplit(rownames(df),"\\."),"[",2)
df$group <- interaction(geno,frac)

# generate the plojt
plot <- ggplot(df, aes(x,y,color=group)) + geom_point(size=4)
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
