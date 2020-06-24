#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: plot grouped protein abundance
#' authors: Tyler W Bradshaw
#' ---

# OPTIONS:

## Input data in root/data/
# * tmt_protein

## Output:
# * a single pdf with all protein plots.
#graphpad_purple <- c("R"=148,"G"=33,"B"=146)
#colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
#    "#942192","#B847B4","#DC6AD7") # Swip Purples

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
ggtheme()
set_font("Arial",font_path=fontdir)

# Load significant modules.
data(sig_modules)

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
## Sort plots by module assignment.
#--------------------------------------------------------------------

# Load the graph partition
data(partition)

# If partition exists, sort the plots by module assignment.
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

#--------------------------------------------------------------------
## Combine plots for all proteins from a module together. 
#--------------------------------------------------------------------

# Module colors.
data(module_colors)

# Module stats for annotations.
data(module_stats)

# Loop to do the work.
grouped_plots <- list()
all_modules <- unique(partition[partition!=0])

for (module in all_modules){
	# Get the protein data for a module.
	prots <- names(which(partition == module))
	module_plots <- plots[prots]
	protein_list <- lapply(module_plots,function(x) x$data)

	# Get modules colors.
	module_color <- module_colors[paste0("M",module)]

	# Define a function that scales things to align.
	norm_to_max <- function(df) {
		df$Normalized.Intensity <- log2(df$Intensity)*(1/max(log2(df$Intensity)))
		return(df)
	}
	# Normalize 
	norm_prot <- lapply(protein_list,norm_to_max)
	prot_df <- bind_rows(norm_prot)

	# Insure Fraction and Cfg forcce are factors.
	# Sort factor levels in a logical order.
	prot_df$Fraction <- factor(prot_df$Fraction,
		      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	prot_df$"Cfg Force (xg)" <- factor(prot_df$"Cfg Force (xg)")
	levels(prot_df$"Cfg Force (xg)") <- c("5,000","9,000","12,000","15,000",
				 "30,000", "79,000","120,000")

	# Fit with lm, add fitted values to df.
	fit <- lm(Normalized.Intensity ~ Fraction + Genotype, data = prot_df)
	prot_df$Fitted.Intensity <- fit$fitted.values

	# Extract other key stats.
	#fit_summary <- summary(fit)
	# Collect some stats.
	#f = paste0("F=",round(fit_summary$fstatistic[1],2))
	#d = paste0("DF=", fit_summary$fstatistic[2])
	#r = paste0("R²=",round(fit_summary$r.squared,3))
	#q = paste0("Adjusted R²=",round(fit_summary$adj.r.squared,3))
	#s = paste0("sigma(σ)= ",round(sigma(fit),3))
	#n <- paste0("n = ",length(unique(prot_df$Accession)))
	#pve <- paste0("PVE = ",
	#	      module_stats %>% 
	#		      filter(Module==as.character(module)) %>% 
	#		      select(PVE) %>% unlist() %>% round(3))
	# Table theme.
	#tab_theme <- ttheme_default()
	#tab_theme$core$fg_params$hjust = 0.5
	#tab_theme$core$bg_params$fill="white"
	#tab_theme$core$bg_params$col=NA
	# Create table.
	#tab <- tableGrob(n, theme=tab_theme, rows=NULL)
	# Add border to table.
	#border <- rectGrob(gp = gpar(fill=NA,lwd=2))
	#gtab <- gtable_add_grob(tab, border, 
	#			t = 1, 
	#			b = nrow(tab), 
	#			l = 1, 
	#			r  = ncol(tab))
	# Generate plot.
	plot <- ggplot(prot_df)
	plot <- plot + aes(x = `Cfg Force (xg)`)
	plot <- plot + aes(y = Normalized.Intensity)
	plot <- plot + aes(group = interaction(Experiment,Treatment,Accession))
	plot <- plot + aes(color = Genotype)
	plot <- plot + geom_line(alpha=0.31)
	plot <- plot + geom_line(aes(y=Fitted.Intensity),size=1.5,alpha=0.5)
	plot <- plot + geom_point(aes(shape=Treatment, fill=Treatment),size=1.5)
	plot <- plot + scale_colour_manual(name="Replicate",
					   values = c(module_color,"#5E5E5E"))
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
	#plot <- plot + 
	#	annotation_custom(gtab, 
	#			  xmin = -Inf, xmax = 2.0, 
	#			  ymin =-Inf, ymax = ymax)

        # Add plot to list.
	grouped_plots[[module]] <- plot
} # EOL
names(grouped_plots) <- paste0("M",c(1:length(grouped_plots)))

# Save.
message("\nSaving all modules.")
myfile <- file.path(root,"figs","Modules","Module_Protein_plots.pdf")
ggsavePDF(grouped_plots, myfile)

# Save significant modules.
message("\nSaving significant modules.")
myfile <- file.path(root,"figs","Modules","Sig_Module_Protein_plots.pdf")
ggsavePDF(grouped_plots[[sig_modules]], myfile)
