#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: plot protein abundance
#' authors: Tyler W Bradshaw
#' ---

# OPTIONS:

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
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

# If necessary, create figsdir.
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# Set theme for the plots:
# Utilize arial font.
set_font("Arial",font_path=fontdir)
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
data(partition)

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

#--------------------------------------------------------------------
## Plot all proteins from a module together. 
#--------------------------------------------------------------------

grouped_plots <- list()

# Loop to do the work.
for (module in seq(1,max(unique(partition)))) {
	# Get  data.
	prots <- names(which(partition == module))
	module_plots <- plots[prots]
	data_list <- lapply(module_plots,function(x) x$data)
	norm_to_max <- function(df) {
		df$Normalized.Intensity <- log2(df$Intensity)*(1/max(log2(df$Intensity)))
		return(df)
	}
	data_list <- lapply(data_list,norm_to_max)
	df <- bind_rows(data_list)
	# Colors for the plot.
	# Generated shade of black online: https://mycolor.space/
	graphpad_purple <- c("R"=148,"G"=33,"B"=146)
	colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
	    "#942192","#B847B4","#DC6AD7") # Swip Purples
	# Insure Fraction is a factor, and levels are in correct order.
	df$Fraction <- factor(df$Fraction,
		      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	df$"Cfg Force (xg)" <- factor(df$"Cfg Force (xg)")
	levels(df$"Cfg Force (xg)") <- c("5,000","9,000","12,000","15,000",
				 "30,000", "79,000","120,000")
	# Fit with glm, add fitted values to df.
	fit <- glm(Normalized.Intensity ~ Fraction + Genotype, data = df)
	df$Fitted.Intensity <- fit$fitted.values
	# Generate plot.
	plot <- ggplot(df)
	plot <- plot + aes(x = Fraction)
	plot <- plot + aes(y = Normalized.Intensity)
	plot <- plot + aes(group = interaction(Experiment,Treatment,Accession))
	plot <- plot + aes(color = Genotype)
	plot <- plot + geom_line(alpha=0.37)
	plot <- plot + geom_line(aes(y=Fitted.Intensity),size=1.5)
	plot <- plot + geom_point(aes(shape=Treatment, fill=Treatment),size=1)
	plot <- plot + scale_colour_manual(name="Replicate",values = c("#DC6AD7","#5E5E5E"))
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
	plot <- plot + theme(axis.text.x = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(axis.text.y = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())
	plot <- plot + theme(legend.position = "none")
	plot <- plot + ggtitle(paste("Module:",module))
	grouped_plots[[module]] <- plot
} # EOL

myfile <- file.path(figsdir,"TMT_Module_Proteins.pdf")
ggsavePDF(grouped_plots, myfile)

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# Save all proteins.
message("\nSaving plots, this will take several minutes.")
myfile <- file.path(figsdir,"TMT_All_Proteins.pdf")
	ggsavePDF(plots, myfile)

message("\nDone!")
