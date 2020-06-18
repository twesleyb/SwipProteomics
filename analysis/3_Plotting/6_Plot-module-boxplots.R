#!/usr/bin/env Rscript

## R Options:

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
myfile <- file.path(root,"rdata","alt_glm_results.RData")
alt_results <- readRDS(myfile)

# Set plotting theme.
ggtheme()
set_font("Arial", font_path = file.path(root,"fonts"))

#---------------------------------------------------------------------
## Prepare the data for plotting.
#---------------------------------------------------------------------

# Annotate tmt data with module membership.
tmt_protein <- tmt_protein %>% filter(Accession %in% names(partition))
tmt_protein$Module <- partition[tmt_protein$Accession]

# Get adjusted fold-changes.
colnames(alt_results)[-c(1,2)] <- paste0("Adjusted.",
					 colnames(alt_results)[-c(1,2)])

# Combine data.
tmt_protein <- left_join(tmt_protein,alt_results,by="Accession")

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
	colors <- c("#47b2a4", module_colors[module_name])
	# Generate the plot.
	plot <- ggplot(df,aes(x=Genotype,y=log2(Intensity),fill=Genotype)) 
	plot <- plot + geom_boxplot() 
	plot <- plot + geom_point() 
	plot <- plot + scale_fill_manual(values=colors)
	plot <- plot + ggtitle(module_name)
	plot <- plot + theme(plot.title = element_text(hjust = 0.5))
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
