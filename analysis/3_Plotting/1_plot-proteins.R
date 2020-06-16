#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: plot protein abundance
#' authors: Tyler W Bradshaw
#' ---

# OPTIONS:
save_subset = TRUE # Only save subset of clustered proteins?

## Input data in root/data/
# * tmt_protein

## Output:
# * a single pdf with all protein plots.

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Load additional functions in root/R/
suppressWarnings({ devtools::load_all() })

# Project directories:
root <- getrd()
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

# If necessary, create figsdir.
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# Utilize arial font.
set_font("Arial",font_path=fontdir)

# Set theme for the plots:
ggtheme()

#--------------------------------------------------------------------
## Generate the plots.
#--------------------------------------------------------------------

# Loop to generate plots for all_proteins.
all_proteins <- unique(tmt_protein$Accession)
message(paste0("\nGenerating plots for all proteins ",
	      "(n=",length(all_proteins),")."))
plots <- list()
pbar <- txtProgressBar(max=length(all_proteins),style=3)
for (prot in all_proteins) {
	plots[[prot]] <- plot_protein(tmt_protein,prot)
	setTxtProgressBar(pbar,value=match(prot,all_proteins))
}
close(pbar)

#--------------------------------------------------------------------
## Save plots as a single pdf.
#--------------------------------------------------------------------

# Load the graph partition
data(list="partition")

# If partition exists, sort the plots by module assignment.
if (exists("partition")) { 
	message("\tSorting plots...")
	# Generate list of modules.
	modules <- split(names(partition),partition)
	names(modules) <- paste("Module:",names(modules))
	# Sort by Module assignment.
	sorted_proteins <- unlist(modules,use.names=FALSE)
	sorted_plots <- plots[sorted_proteins]
	remainder_prots <- names(plots)[names(plots) %notin% names(sorted_plots)]
	remainder_plots <- plots[remainder_prots]
	plots <- c(sorted_plots,remainder_plots)
	# Annotate plots with module assignment.
	message("\tAnnotating plots...")
	for (i in c(1:length(plots))) {
		protein <- names(plots)[i]
		plot <- plots[[protein]]
		module <- paste("Module:", partition[protein])
		yrange <- plot$data %>% filter(Accession == protein) %>% 
			select(Intensity) %>% log2() %>% range()
		ypos <- yrange[1] - 0.1* diff(yrange)
		plot <- plot + annotate(geom="label",x=7, y=ypos, label=module)
		plots[[protein]] <- plot
	}
} # EIF

# Save.
message("\nSaving plots, this will take several minutes.")
myfile <- file.path(figsdir,"TMT_All_Proteins.pdf")
if (save_subset) {
	# Save subset of clustered proteins.
	idx <- names(plots)[names(plots) %in% names(partition)]
	ggsavePDF(plots[idx], myfile)
} else {
	# Save all proteins.
	ggsavePDF(plots, myfile)
}

message("\nDone!")
