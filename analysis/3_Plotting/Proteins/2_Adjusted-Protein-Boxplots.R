#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## INPUTS:
#* tmt_data

## OPTIONS:
fig_width = 5
fig_height =  5
colors = c(WT="#47b2a4",MUT="#B86FAD") # Colors for WT and mutant groups.

## OUTPUTS:
# Boxplots of protein abundance adjusted for baseline differences in
# fraction.

#--------------------------------------------------------------------
## Misc function - getrd
#--------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
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

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE)

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
figsdir <- file.path(rootdir, "figs","Proteins") 
suppdir <- file.path(rootdir,"manuscript","files")

# Set global plotting settings.
ggtheme(); set_font("Arial", font_path = fontdir)

# Load the data.
data(gene_map)
data(tmt_protein)

#---------------------------------------------------------------------
## Plots for select proteins.
#---------------------------------------------------------------------

# Wash proteins + control.
prots <- c("Washc4","Washc1","Washc2","Washc5","Tubb4a")
idx <- match(prots,gene_map$symbol)
names(prots) <- gene_map$uniprot[idx]

# Subset the data.
df <- tmt_protein %>% filter(Accession %in% names(prots))

# Labels will simply be WT and Mutant.
xlabels <- rep(c("WT","MUT"),times=length(prots))

# Order of the factors.
df$Group <- df$Treatment
df$Group <- gsub("Control","WT",df$Group)
df$Group <- gsub("Mutant","MUT",df$Group)
factor_order <- paste(rep(names(prots),each=2), c("WT","MUT"),sep=".")
df$Accession.Group <- as.character(interaction(df$Accession,df$Group))
df$Accession.Group <- factor(df$Accession.Group,levels=factor_order)

# Generate a plot.
plot <- ggplot(df, aes(x=Accession.Group, y=log2(Adjusted.Intensity),
		       fill=Group)) + 
	geom_boxplot() + geom_point(aes(fill=Group,shape=Group)) +
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
message("\nSelect Protein GLM stats:")
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
build <- ggplot_build(plot)
ymax <- build$layout$panel_params[[1]][["y.range"]][2]
plot <- plot + annotate("text",x=seq(1.5,length(prots)*2,by=2),
			y=ymax,label=symbols,size=5)

# Save as pdf.
myfile <- file.path(figsdir,"WASH_Complex_Protein_BoxPlots.pdf")
ggsave(myfile,plot, height = fig_height, width = fig_width)

#----------------------------------------------------------------------
## Plot all sig proteins.
#----------------------------------------------------------------------

# Any protein with FDR < 0.1
data(sig_proteins)

# Loop to do the work.
plot_list <- list()
for (prot in sig_proteins) {
	# Subset the data.
	df <- tmt_protein %>% filter(Accession == prot)
	# Organize the factors.
	df$Group <- df$Genotype
	df$Group <- factor(df$Group, levels= c("WT","MUT"))
	# Generate a plot.
	plot <- ggplot(df, aes(x=Group, y=log2(Adjusted.Intensity),fill=Group)) 
	plot <- plot + geom_boxplot() 
	plot <- plot + geom_point(aes(fill=Group,shape=Group))
	#plot <- plot + theme(axis.text.x=element_text(angle=45))
	plot <- plot + scale_fill_manual(name="Genotype",values=colors)
	plot <- plot + theme(legend.position="none")
	plot <- plot + theme(panel.background = element_blank())
	#plot <- plot + theme(panel.border=element_rect(colour="black",fill="NA",size=1)) # border around entire plot
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())
	plot <- plot + theme(axis.title.x = element_blank())
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
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
	idx <- match(prot,gene_map$uniprot)
	symbol <- gene_map$symbol[idx]
	plot <- plot + ggtitle(paste0(symbol,"|",prot))
	# Add module annotation.
	yrange <- log2(range(df$Adjusted.Intensity))
	ypos <- max(yrange) + 0.05 * diff(yrange)
	mylabel <- paste("Module:",partition[prot])
	plot <- plot + annotate(geom="label",x=2.25, y=ypos, label=mylabel)
	plot_list[[prot]] <- plot
}

# Save as pdf.
message("\nSaving protein boxplots.")
myfile <- file.path(suppdir,"S5_Adjusted_Protein_BoxPlots.pdf")
ggsavePDF(plot_list,myfile)
