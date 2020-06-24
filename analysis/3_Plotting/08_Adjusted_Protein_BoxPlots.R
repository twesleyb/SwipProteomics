#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## INPUTS:

## OPTIONS:
fig_width = 5
fig_height =  5
colors = c(WT="#47b2a4",MUT="#B86FAD") # Colors for WT and mutant groups.

## OUTPUTS:
# Boxplot of select protein abundance adjusted for baseline differences in
# fraction.
# Boxplots of all proteins.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE) # NOTE: getrd is a f(x) in .Rprofile.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(ggplot2) # For making plots.
	library(data.table) # For working with tables.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories:
fontdir <- file.path(rootdir, "fonts") # Arial font for plots.
figsdir <- file.path(rootdir, "figs","Proteins") # Output figures.

# Set global plotting settings.
ggtheme()
set_font("Arial", font_path = fontdir)

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

data(gene_map)
data(tmt_protein)

# Wash proteins + control.
prots <- c("Washc4","Washc1","Washc2","Washc5","Tubb4a")
names(prots) <- gene_map$"Uniprot Accession"[match(prots,gene_map$"Gene Symbol")]

# Subset the data.
df <- tmt_protein %>% filter(Accession %in% names(prots))

# Labels will simply be WT and Mutant.
xlabels <- rep(c("WT","Mutant"),times=length(prots))

# Order of the factors.
factor_order <- paste(rep(names(prots),each=2), c("Control","Mutant"),sep=".")
df$Accession.Treatment <- as.character(interaction(df$Accession,df$Treatment))
df$Accession.Treatment <- factor(df$Accession.Treatment,levels=factor_order)

# Generate a plot.
plot <- ggplot(df, aes(x=Accession.Treatment, y=log2(Adjusted.Intensity),
		       fill=Treatment)) + 
	geom_boxplot() + geom_point(aes(fill=Treatment,shape=Treatment)) +
	theme(axis.text.x=element_text(angle=45))
plot <- plot + scale_fill_manual(name="Genotype",values=colors)
plot <- plot + theme(legend.position="none")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border=element_rect(colour="black",fill="NA",size=1))
plot <- plot + theme(axis.title.x = element_blank())
plot <- plot + scale_x_discrete(labels=xlabels)
plot <- plot + ylab("log2(Adjusted Intensity)")

# Add some lines to break up the data.
plot <- plot + geom_vline(xintercept=seq(2.5,length(prots)*2,by=2),
			  linetype="dotted",size=0.5)

# Pretty print select protein stats:
message("\nSelect Protein stats:")
stats <- df %>% filter(Genotype == "MUT") %>%
	select(Accession, Symbol, Adjusted.logFC, Adjusted.PercentWT, 
	       Adjusted.F, Adjusted.FDR) %>% unique()
knitr::kable(stats)

# Significance stars:
stats$symbol <- "ns"
stats$symbol[stats$Adjusted.FDR < 0.05] <- "*"
stats$symbol[stats$Adjusted.FDR < 0.005] <- "**"
stats$symbol[stats$Adjusted.FDR < 0.0005] <- "***"
stats$xpos <- seq(1.5,by=2,length.out=5)
yrange <- range(log2(df$Adjusted.Intensity))
stats$ypos <- max(yrange) + 0.05 * diff(yrange)

# Add significance stars.
plot <- plot + 
	annotate("text",x=stats$xpos,y=stats$ypos,label=stats$symbol,size=4)

# Annotate with protein names.
symbols <- prots 
#----------------------------------------------------------------------
build <- ggplot_build(plot)
ymax <- build$layout$panel_params[[1]][["y.range"]][2]
plot <- plot + annotate("text",x=seq(1.5,length(prots)*2,by=2),
			y=ymax,label=symbols,size=5)

# Save as pdf.
myfile <- file.path(figsdir,"Select_Adjusted_Protein_BoxPlots.pdf")
ggsave(myfile,plot, height = fig_height, width = fig_width)

#----------------------------------------------------------------------
## Plot all sig proteins.
#----------------------------------------------------------------------

# Any protein with FDR < 0.1
data(sig_proteins)
prots <- sig_proteins$sig968

# Loop to do the work.
plot_list <- list()
for (protein in prots) {
	# Subset the data.
	df <- tmt_protein %>% filter(Accession == protein)
	# Organize the factors.
	df$Group <- df$Treatment
	df$Group <- gsub("Control","WT",df$Group)
	df$Group <- gsub("Mutant","MUT",df$Group)
	df$Group <- factor(df$Group)
	levels(df$Group) <- c("WT","MUT")
	# Generate a plot.
	plot <- ggplot(df, aes(x=Group, y=log2(Adjusted.Intensity),fill=Group)) 
	plot <- plot + geom_boxplot() 
	plot <- plot + geom_point(aes(fill=Group,shape=Group))
	plot <- plot + theme(axis.text.x=element_text(angle=45))
	plot <- plot + scale_fill_manual(name="Genotype",values=colors)
	plot <- plot + theme(legend.position="none")
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(panel.border=element_rect(colour="black",
						       fill="NA",size=1))
	plot <- plot + theme(axis.title.x = element_blank())
	plot <- plot + scale_x_discrete(labels=xlabels)
	plot <- plot + ylab("log2(Adjusted Intensity)")
	# Collect protein stats.
	stats <- df %>% filter(Group == "MUT") %>%
		select(Accession, Symbol, Adjusted.logFC, Adjusted.PercentWT, 
		       Adjusted.F, Adjusted.FDR) %>% unique()
	# Significance stars:
	stats$symbol <- "ns"
	stats$symbol[stats$Adjusted.FDR < 0.05] <- "*"
	stats$symbol[stats$Adjusted.FDR < 0.005] <- "**"
	stats$symbol[stats$Adjusted.FDR < 0.0005] <- "***"
	stats$xpos <- 1.5
	yrange <- range(log2(df$Adjusted.Intensity))
	stats$ypos <- max(yrange) + 0.05 * diff(yrange)
	# Add significance stars.
	plot <- plot + 
		annotate("text",x=stats$xpos,y=stats$ypos,
			 label=stats$symbol,size=4)
	# Add title.
	idx <- match(protein,gene_map$"Uniprot Accession")
	symbol <- gene_map$"Gene Symbol"[idx]
	plot <- plot + ggtitle(paste0(symbol,"|",protein))
	# Add module annotation.
	yrange <- log2(range(df$Adjusted.Intensity))
	ypos <- max(yrange) + 0.05 * diff(yrange)
	mylabel <- paste("Module:",partition[protein])
	plot <- plot + annotate(geom="label",x=2.25, y=ypos, label=mylabel)
	plot_list[[protein]] <- plot
}

# Save as pdf.
message("\nSaving protein boxplots.")
myfile <- file.path(figsdir,"Sig968_Adjusted_Protein_BoxPlots.pdf")
ggsavePDF(plot_list,myfile)
