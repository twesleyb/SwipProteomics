#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: plot protein abundance
#' authors: Tyler W Bradshaw
#' ---

## Options:
save_all = FALSE
save_sig = TRUE
FDR_alpha = 0.1

## Input data in root/data/
# * tmt_protein

## Output:
# * a single pdf with plots of all proteins

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
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

# If necessary, create figsdir.
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# Set theme for the plots:
ggtheme(); set_font("Arial",font_path=fontdir)

# Load sig85 proteins.
data(tmt_protein)
sig_prots <- tmt_protein %>% filter(FDR < FDR_alpha) %>% 
	select(Accession) %>% unlist() %>% unique()

#--------------------------------------------------------------------
## Generate the plots.
#--------------------------------------------------------------------
# Loop to generate plots for all_proteins.

plots <- list()
all_proteins <- unique(tmt_protein$Accession)
message(paste("\nGenerating plots for",
	       formatC(length(all_proteins),big.mark=","),
	       "proteins."))
pbar <- txtProgressBar(max=length(all_proteins),style=3)
for (prot in all_proteins) {
	plots[[prot]] <- plot_protein(tmt_protein,prot)
	setTxtProgressBar(pbar,value=match(prot,all_proteins))
}
close(pbar)

# Generate a plot with a legend.
plot <- plot_protein(tmt_protein,"P08226",legend=TRUE)
plot_legend <- cowplot::get_legend(plot) 
plot  <- plot + theme(legend.position="none")
plot <- plot + theme(axis.text.x = element_text(angle = 45))
plot <- cowplot::plot_grid(plot,plot_legend,rel_widths=c(1,0.2))
myfile <- file.path(root,"figs","Examples","S4_Example.png")
ggsave(plot,file=myfile,width=4.5,height=4.5)

#--------------------------------------------------------------------
## Sort plots by module membership and save as a single pdf.
#--------------------------------------------------------------------
# NOTE: Proteins that are not clustered (M0), are not saved.

# Load the graph partition.
data(partition)

# Generate list of modules.
modules <- split(names(partition),partition)
names(modules) <- paste("Module:",names(modules))

# Sort plots by Module assignment.
sorted_proteins <- unlist(modules,use.names=FALSE)
sorted_plots <- plots[sorted_proteins]
plots <- sorted_plots

# Annotate plots with module assignment.
for (i in c(1:length(plots))) {
	protein <- names(plots)[i]
	plot <- plots[[protein]]
	module <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Accession == protein) %>% 
		select(Intensity) %>% log2() %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plot <- plot + annotate(geom="label",x=7, y=ypos, label=module)
	plots[[protein]] <- plot
}

#--------------------------------------------------------------------
## Save a single plot as an example.
#--------------------------------------------------------------------

# Example protein.
set.seed(7)
prot <- sample(sig_prots,1)
plot <- plots[[prot]]
myfile <- file.path(figsdir,"S4_Example.png")
ggsave(plot,file=myfile,height=7,width=7)


#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

if (save_all) {
	# Save all proteins.
	message("\nSaving all plots, this will take several minutes.")
	M0_prots <- names(partition[partition==0])
	drop <- names(plots) %in% M0_prots
	all_plots <- plots[!drop]
	myfile <- file.path(figsdir,"Protein_plots.pdf")
	ggsavePDF(all_plots, myfile)
}

# Save sig85 prots include proteins that were assigned to M0.
if (save_sig) {
	message("\nSaving significant plots, this will take several minutes.")
	myfile <- file.path(figsdir,"Sig85_Protein_plots.pdf")
	prots <- sig_proteins[["sig85"]]
	# Sort by module assignment.
	idx <- order(partition[prots])
	sig_plots <- plots[prots][idx]
	ggsavePDF(sig_plots, myfile)
}
