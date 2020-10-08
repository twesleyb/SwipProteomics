#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
#save_all = FALSE
#save_sig = TRUE
FDR_alpha = 0.1
PAdj_alpha = 0.05

# Input data in root/data/
# * tmt_protein


## Output ---------------------------------------------------------------------
# * a single pdf with plots of all proteins


## Misc function - getrd ------------------------------------------------------

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

## Set-up the workspace.

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

# load the data
data(msstats_prot)

############################################################################
# WORK

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

	# Colors for the plot.
	# Generated shade of black online: https://mycolor.space/
	graphpad_purple <- c("R"=148,"G"=33,"B"=146)
	colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
	   	    "#942192","#B847B4","#DC6AD7") # Swip Purples
	
	# Subset the data.
	df <- subset(data,data[[id.col]] == protein)
	gene <- unique(df$Symbol)

	# Insure Fraction is a factor, and levels are in correct order.
	df$Fraction <- factor(df$Fraction,
			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	df$"Cfg Force (xg)" <- factor(df$"Cfg Force (xg)")
	levels(df$"Cfg Force (xg)") <- c("5,000","9,000","12,000","15,000",
					 "30,000", "79,000","120,000")

	# Collect FDR stats.
	stats <- df %>% group_by(Treatment,Fraction) %>% 
		summarize(Intensity = max(Intensity), FDR = unique(FDR))
	stats$ypos <- 1.02 * log2(max(stats$Intensity))
	stats <- stats %>% filter(Treatment == "Control")
	stats$Cfg.Force <- levels(df$"Cfg Force (xg)")
	stats$symbol <- ""
	stats$symbol[stats$FDR<0.1] <- "."
	stats$symbol[stats$FDR<0.05] <- "*"
	stats$symbol[stats$FDR<0.005] <- "**"
	stats$symbol[stats$FDR<0.0005] <- "***"
	#stats$ypos <- 1.02*log2(stats$Intensity)

	# Generate the plot.
	plot <- ggplot(df, aes(x = Fraction, y = log2(Intensity),
			       group = interaction(Experiment,Treatment),
			       colour = interaction(Experiment,Treatment))) + 
                geom_point(aes(shape=Treatment,
			       fill=Treatment),size=2) + 
		geom_line() + 
		ggtitle(protein)
	# Annotate with significance stars.
	plot <- plot + annotate("text",x=stats$Fraction,
				y=max(stats$ypos),label=stats$symbol,size=7)
	# Add Custom colors and modify legend title and labels.
	mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
	plot <- plot + scale_colour_manual(name="Replicate",
				           values=colors,
					   labels=mylabs) 
	plot <- plot + scale_x_discrete(labels=stats$Cfg.Force)
	plot <- plot + xlab("Force (xg)")
	plot <- plot + ggtitle(paste(gene,protein,sep=" | "))

	# Edit y axis.
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))

	# Make x and y axes consistent.
	plot <- plot + theme(axis.text.x = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(axis.text.y = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))

	# Remove background.
	plot <- plot + theme(panel.background = element_blank())

	# Add x and y axes.
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())

	# Remove legend.
	if (!legend) {
		plot <- plot + theme(legend.position = "none")
	}

	return(plot)
}



############################################################################
## Generate the plots.

# Loop to generate plots for all_proteins.

proteins <- unique(as.character(msstats_prot$Protein))

plots <- list()

message("\nGenerating plots for ", 
	formatC(length(proteins),big.mark=","), " proteins.")
pbar <- txtProgressBar(max=length(proteins),style=3)

for (prot in proteins) {

	plots[[prot]] <- plot_protein(ms,prot)

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
