#!/usr/bin/env Rscript

# Load renv.
root <- getrd()
renv::load(root)

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Project imports.
devtools::load_all()

# Load the proteomics data.
data(tmt_protein)

# Load the graph partition.
data(partition)

# Load module color assignments.
data(module_colors)

# Load statistics from adjusted data.
data(tmt_protein)

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
modules <- modules[order(modules)][-1]

# Loop to generate plots.
plots <- list()
for (module in modules) {

	# Get data for a given module.
	module_name <- paste0("M",module)
	df <- tmt_protein %>% filter(Module == module) %>% 
		group_by(Module,Genotype,Fraction) %>%
		summarize(Intensity=sum(2^Adjusted.Intensity),.groups="drop")
	df$Genotype <- factor(df$Genotype,levels=c("WT","MUT"))

	# Get module color.
	colors <- c(col2hex("gray"), module_colors[module_name])

	# Collect the module stats.
	stats <- module_stats %>% filter(Module == module) %>% 
		select(Module, Nodes, PVE, Hubs)
	
	## FIXME: annotate with Module stats.

    yrange <- unlist(build$layout$panel_params[[1]][8])
    xrange <- range(plot$data$Log2QC1)
    xmin <- min(xrange)
    xmax <- max(xrange)
    xdelta <- xmax - xmin
    ymin <- min(yrange)
    ymax <- max(yrange)
    ydelta <- ymax - ymin
    tt <- ttheme_default(
      base_size = 11,
      core = list(bg_params = list(fill = "white"))
    )
    tab <- tableGrob(mytable, rows = NULL, theme = tt)
    g <- gtable_add_grob(tab,
      grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
      t = 1, b = nrow(tab), l = 1, r = ncol(tab)
    )

    if (annotate == TRUE) {
      plot <- plot + annotation_custom(g,
        xmin = xmin - 0.65 * xdelta, xmax,
        ymin = ymin + 0.8 * ydelta, ymax
      )
    }

	# Generate the plot.
	plot <- ggplot(df,aes(x=Genotype,y=log2(Intensity),fill=Genotype)) 
	plot <- plot + geom_boxplot() 
	plot <- plot + geom_point() 
	plot <- plot + scale_fill_manual(values=colors)
	plot <- plot + ggtitle(module_name)
	plot <- plot + theme(plot.title = element_text(hjust = 0.5))
	plot <- plot + theme(legend.position = "none")

	y <- max(log2(plot$data$Intensity))
	plot + geom_label(x=1.5,y=y,label="foo",fill="NA")


	plot + ggtitle(paste(module_name,"\nPVE:",round(stats$PVE,3)))

# N Nodes, PVE, HubS
	plot + geom_label(x=1.5,y=y,label="",fill="white")




	
	module_stats$Module


	plots[[module_name]] <- plot
}

# FIXME: annotate with pvalue.
# FIXME: Insure wash module is correct color. (Mutant="#B86FAD")
# FIXME: change plot points -- mutant is triangle
# FIXME: clean up plot axes and background

# Save as single pdf.
# NOTE: This takes a couple minutes.
myfile <- file.path(root,"figs","Modules","All_Modules.pdf")
ggsavePDF(plots, myfile)
