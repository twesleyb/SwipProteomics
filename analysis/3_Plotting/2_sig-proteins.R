#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: plot protein abundance
#' authors: Tyler W Bradshaw
#' ---

## Input data in root/data/
# * tmt_protein

## Output:
# * a single pdf with plots of significantly DA proteins.

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
if (!dir.exists(figsdir))  {
	dir.create(figsdir)
}

# Utilize arial font.
set_font("Arial",font_path=fontdir)

# Set theme for the plots:
ggtheme()

#--------------------------------------------------------------------
## Prepare the data for plotting.
#--------------------------------------------------------------------

# Load the normalized protein data.
data(tmt_protein)

# Collect proteins with any sig change.
sig_prots <- tmt_protein %>% filter(FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()

#--------------------------------------------------------------------
## Generate the plots.
#--------------------------------------------------------------------

# Loop to generate all plots.
message(paste0("\nGenerating plots of significant proteins ",
	      "(n=",length(sig_prots),")."))
plots <- list()
pbar <- txtProgressBar(max=length(sig_prots),style=3)
for (protein in sig_prots) {
	plots[[protein]] <- plot_protein(tmt_protein,protein)
	setTxtProgressBar(pbar,value=match(protein,sig_prots))
}
close(pbar)

#--------------------------------------------------------------------
## Save all significant proteins as a single pdf.
#--------------------------------------------------------------------

# Load the graph partition
data(list=partition) # FIXME: warnings about datasets?

# If partition exists, sort the plots by module assignment.
if (exists("partition")) { 
	message("\tSorting plots...")
	# Generate list of modules.
	modules <- split(names(partition),partition)
	names(modules) <- paste("Module:",names(modules))
	# Sort by Module assignment.
	sorted_proteins <- unlist(modules,use.names=FALSE)
	sorted_proteins <- sorted_proteins[sorted_proteins %in% names(plots)]
	sorted_plots <- plots[sorted_proteins]
	# Annotate plots with module assignment.
	message("\tAnnotating plots...")
	for (i in c(1:length(sorted_plots))) {
		protein <- names(sorted_plots)[i]
		plot <- sorted_plots[[protein]]
		module <- paste("Module:", partition[protein])
		yrange <- plot$data %>% filter(Accession == protein) %>% 
			select(Intensity) %>% log2() %>% range()
		ypos <- yrange[1] - 0.1* diff(yrange)
		plot <- plot + annotate(geom="label",x=7, y=ypos, label=module)
		sorted_plots[[protein]] <- plot
	}
} # EIF

# Save significant proteins as a single pdf.
myfile <- file.path(figsdir,"TMT_Significant_Proteins.pdf")
ggsavePDF(sorted_plots,myfile)

message("\nDone!")
