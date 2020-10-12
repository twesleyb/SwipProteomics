#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: generate module plots showing all proteins
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:
save_all = TRUE
save_sig = FALSE
wt_color = "#47b2a4" # teal blue

## Input data in root/data/
# * tmt_protein

## Output:
# * a single pdf with the aligned protein plots for all modules. 

#---------------------------------------------------------------------
## Misc functions
#---------------------------------------------------------------------

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

# Define a function that scales things to align all proteins from a
# module.
norm_to_max <- function(df) {
	norm <- function(x) { log2(x)*(1/max(log2(x))) }
	df$Normalized.Intensity <- norm(df$Intensity)
	return(df)
}

# Get median replicate for each protein.
# Using the median simplifies the appearance of thee plots.
get_median <- function(x) {
	df <- x %>% group_by(Experiment,Fraction,
			     Treatment,Accession) %>% 
		summarize("Median(Intensity)"=median(Intensity),
			  .groups="drop")
	return(df)
}

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	# For working with data tables:
	library(dplyr)
	library(data.table)
	# For plotting:
	library(ggplot2)
	# For working with tables as graphics:
	library(grid)
	library(gtable)
	library(gridExtra)
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

# Set theme for the plots; utilize arial font.
ggtheme(); set_font("Arial",font_path=fontdir)

# Load significant modules.
data(sig_modules) # 37 sig modules

#--------------------------------------------------------------------
## Generate the plots.
#--------------------------------------------------------------------

# Problem protein: "E9Q6J5" 
#partition["E9Q6J5"] = M222

# All proteins. Remove problematic protein.
all_proteins <- unique(tmt_protein$Accession)
#all_proteins <- all_proteins[all_proteins!="E9Q6J5"]

# Loop to generate plots for all_proteins.
message(paste0("\nGenerating plots for all proteins..."))
plots <- list()
pbar <- txtProgressBar(max=length(all_proteins),style=3)
for (prot in all_proteins) {
	plots[[prot]] <- plot_protein(tmt_protein,prot)
	setTxtProgressBar(pbar,value=match(prot,all_proteins))
}
close(pbar)

#--------------------------------------------------------------------
## Sort plots by module assignment.
#--------------------------------------------------------------------

# Load the graph partition
data(partition)

# If partition exists, sort the plots by module assignment.
message("\tSorting plots by module assignment.")

# Generate list of modules.
modules <- split(names(partition),partition)
names(modules) <- paste("Module:",names(modules))

# Sort by Module assignment.
sorted_proteins <- unlist(modules,use.names=FALSE)
sorted_plots <- plots[sorted_proteins]
remainder_prots <- names(plots)[names(plots) %notin% names(sorted_plots)]
remainder_plots <- plots[remainder_prots]
plots <- c(sorted_plots,remainder_plots)

# Drop pesky prot.
nout <- sum(is.na(names(plots)))
message(paste("Warning: removing",nout,"plot(s)."))
plots <- plots[-which(is.na(names(plots)))]

# Annotate plots with module assignment.
message("\tAnnotating plots with module assignment.")
pbar <- txtProgressBar(max=length(plots),style=3)
for (i in c(1:length(plots))) {
	protein <- names(plots)[i]
	plot <- plots[[protein]]
	module <- paste("Module:", partition[protein])
	yrange <- plot$data %>% filter(Accession == protein) %>% 
		select(Intensity) %>% log2() %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plot <- plot + annotate(geom="label",x=7, y=ypos, label=module)
	plots[[protein]] <- plot
	setTxtProgressBar(pbar,value=i)
}
close(pbar)

#--------------------------------------------------------------------
## Combine plots for all proteins from a module together. 
#--------------------------------------------------------------------

# Module colors.
data(module_colors)

# Loop to do the work.
grouped_plots <- list()
all_modules <- unique(partition[partition!=0])
all_modules <- all_modules[order(all_modules)] # Sort

# Loop to do work.
message("\nAligning module prortein plots.")
pbar <- txtProgressBar(max=length(all_modules),style=3)
for (module in all_modules){
	# Get the protein data for a module.
	prots <- names(which(partition == module))
	module_plots <- plots[prots]
	protein_list <- lapply(module_plots,function(x) x$data)
	# Remove any null data arising from pesky prot.
	idx <- sapply(protein_list,is.null)
	if (any(idx)) { 
		message(paste("Warning: NULL data is removed."))
		protein_list <- protein_list[-which(idx)] 
	}
	# Get modules colors.
	module_color <- module_colors[paste0("M",module)]
	# Normalize/scale.
	norm_prot <- lapply(protein_list,norm_to_max)
	# Combine data for all proteins together.
	prot_df <- bind_rows(norm_prot)
	# To simplify plot, calculate protein-wise mean.
	prot_df <- prot_df %>% group_by(Accession,`Cfg Force (xg)`,
					Fraction,Treatment) %>%
		summarize(Normalized.Intensity=mean(Normalized.Intensity),
			  .groups="drop")
	# Insure Fraction and Cfg force are factors.
	# Sort factor levels in a logical order.
	prot_df$Fraction <- factor(prot_df$Fraction,
		      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	prot_df$"Cfg Force (xg)" <- factor(prot_df$"Cfg Force (xg)")
	levels(prot_df$"Cfg Force (xg)") <- c("5,000","9,000","12,000","15,000",
				 "30,000", "79,000","120,000")
	# Fit with lm, add fitted values to df.
	fit <- lm(Normalized.Intensity ~ Fraction + Treatment, data = prot_df)
	prot_df$Fitted.Intensity <- fit$fitted.values
	# Generate plot.
	plot <- ggplot(prot_df)
	plot <- plot + aes(x = `Cfg Force (xg)`)
	plot <- plot + aes(y = Normalized.Intensity)
	plot <- plot + aes(group = interaction(Treatment,Accession))
	plot <- plot + aes(color = Treatment)
	plot <- plot + geom_line(alpha=0.31)
	plot <- plot + geom_line(aes(y=Fitted.Intensity),size=1.5,alpha=0.5)
	plot <- plot + geom_point(aes(shape=Treatment, fill=Treatment),size=1.5)
	plot <- plot + scale_colour_manual(name="Replicate",
					   values = c(wt_color,module_color))
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
	# Add module annotations.
	yrange <- range(plot$data$Normalized.Intensity)
	ymax <- yrange[1] + 0.10 * diff(yrange)
        # Add plot to list.
	grouped_plots[[module]] <- plot
	setTxtProgressBar(pbar,value=match(module,all_modules))
} # EOL
close(pbar)
names(grouped_plots) <- paste0("M",c(1:length(grouped_plots)))

# Save.
if (save_all) {
	message("\nSaving all modules, this will take several minutes.")
	myfile <- file.path(root,"figs","Modules","Module_Protein_plots.pdf")
	ggsavePDF(grouped_plots, myfile)
}

# Save significant modules.
if (save_sig) {
	message("\nSaving significant modules.")
	myfile <- file.path(figsdir, 
			    paste0("Sig",length(sig_modules),"_Module_Protein_plots.pdf"))
	ggsavePDF(grouped_plots[sig_modules], myfile)
}
