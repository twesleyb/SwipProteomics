#!/usr/bin/env Rscript

## OPTIONS:
BF_alpha = 0.05

# Load renv.
root <- getrd()
renv::load(root)

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
	# For working with tables as graphics:
	library(grid)
	library(gtable)
	library(gridExtra)
})

# Project imports.
devtools::load_all()

# Load the proteomics data.
data(tmt_protein)

# Load the graph partition.
data(partition)

# Load module color assignments.
data(module_colors)

# Load the module-level statistics.
data(module_stats)

# Set plotting theme.
ggtheme()
set_font("Arial", font_path = file.path(root,"fonts"))

#---------------------------------------------------------------------
## Prepare the data for plotting.
#---------------------------------------------------------------------

# Annotate tmt data with module membership.
tmt_protein$Module <- partition[tmt_protein$Accession]

# All modules.
modules <- unique(tmt_protein$Module)
modules <- modules[order(modules)][-1] # Sort in numerical order.

# Loop to generate plots.
plots <- list()
for (module in modules) {

	# Get data for a given module.
	module_name <- paste0("M",module)
	df <- tmt_protein %>% filter(Module == module) %>% 
		group_by(Module,Genotype,Fraction) %>%
		summarize(Intensity=sum(Adjusted.Intensity),.groups="drop")
	df$Genotype <- factor(df$Genotype,levels=c("WT","MUT"))

	# Get module's color.
	colors <- c(col2hex("gray"), module_colors[module_name])

	# Collect the module's stats.
	stats <- module_stats %>% filter(Module == module) %>% 
		select(Module, Nodes, PVE, Hubs, PAdjust)

	# Significance annotations.
	stats$symbol <- ""
	if (stats$PAdjust < 0.1) { stats$symbol <- "." }
	if (stats$PAdjust < 0.05) { stats$symbol <- "*" }
	if (stats$PAdjust < 0.005) { stats$symbol <- "**" }
	if (stats$PAdjust < 0.0005) { stats$symbol <- "***" }

	# Generate the plot.
	plot <- ggplot(df,aes(x=Genotype,y=log2(Intensity),fill=Genotype)) 
	plot <- plot + geom_boxplot() 
	plot <- plot + geom_point(aes(shape=Genotype, fill=Genotype),size=2)
	plot <- plot + scale_fill_manual(values=colors)
	plot <- plot + ggtitle(module_name)
	plot <- plot + theme(plot.title = element_text(hjust = 0.5))
	plot <- plot + theme(legend.position = "none")
	plot <- plot + ylab("Log2(Adjusted Module Intensity)")
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))

	# Add significance star.
	if (stats$PAdjust < BF_alpha) {
		yrange <- range(log2(df$Intensity))
		ypos <- yrange[2] + 0.05 * diff(yrange)
		plot <- plot + 
			annotate("text",x=1.5,size=7,
				 y=ypos, label=stats$symbol)
	}
		
	## Add a table with module statistics.
	# Table theme.
	tab_theme <- ttheme_default()
	tab_theme$core$fg_params$hjust = 0.5
	tab_theme$core$bg_params$fill="white"
	tab_theme$core$bg_params$col=NA
	
	# Create table.
	idx <-  module_stats$Module == as.character(module)
	pve <- paste0("PVE=",round(module_stats$PVE[idx],3))
	padj <- paste0("P-adjust=",round(stats$PAdjust,3))
	n <- paste0("n Proteins=",module_stats$Nodes[idx])
	tab <- tableGrob(rbind(padj,pve,n), theme=tab_theme, rows=NULL)

	# Add border to table.
	border <- rectGrob(gp = gpar(fill=NA,lwd=2))
	gtab <- gtable_add_grob(tab, border, 
				t = 1, 
				b = nrow(tab), 
				l = 1, 
				r  = ncol(tab))

	# Add table to plot.
	# Should table be positioned on the left or right?
	min_group <- df %>% ungroup() %>% 
		filter(Intensity == min(Intensity)) %>% 
		select(Genotype) %>% unlist() %>%
		as.character()
	if (min_group == "WT") {
		xpos <- 4
	} else {
		xpos <- 1
	}
	yrange <- range(log2(df$Intensity))
	ypos <- yrange[1] + 0.15* diff(yrange)
	plot <- plot + 
		annotation_custom(gtab, xmin=-Inf,xmax=xpos,ymin=-Inf,ymax=ypos)

	plots[[module_name]] <- plot
}

# Save as single pdf.
# NOTE: This takes a couple minutes.
myfile <- file.path(root,"figs","Modules","Module_Boxplots.pdf")
ggsavePDF(plots, myfile)
