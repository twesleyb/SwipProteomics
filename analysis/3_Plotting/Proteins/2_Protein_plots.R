#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
FDR_alpha = 0.05
PAdj_alpha = 0.05

# Input data in root/data/
# * msstats_prot


## Output ---------------------------------------------------------------------
# * a single pdf with plots of all proteins


## Functions -------------------------------------------------------------------

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


plot_protein <- function(prot_df, gene_map, protein,legend=FALSE) {
  # function to generate a proteins plot
  # Subset the data.
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  df <- subset(prot_df,prot_df$Protein == protein)
  # Insure Fraction is a factor, and levels are in correct order.
  df$BioFraction <- factor(df$BioFraction,
  			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Collect FDR stats.
  stats <- df %>% group_by(Genotype,BioFraction) %>% 
  		summarize(`Max Abundance` = max(Abundance),
  			  FDR = unique(adj.pvalue),.groups="drop")
  stats$ypos <- 1.02 * max(stats$`Max Abundance`)
  stats <- stats %>% filter(Genotype == "Control")
  stats$symbol <- ""
  stats$symbol[stats$FDR<0.1] <- "."
  stats$symbol[stats$FDR<0.05] <- "*"
  stats$symbol[stats$FDR<0.005] <- "**"
  stats$symbol[stats$FDR<0.0005] <- "***"
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction, y = Abundance, 
  		   group = interaction(Mixture,Genotype),
  		   colour = interaction(Mixture,Genotype))
  plot <- plot + geom_point(aes(shape=Genotype, fill=Genotype),size=2)
  plot <- plot + geom_line()
  plot <- plot + ggtitle(protein)
  # Annotate with significance stars.
  if (any(stats$FDR<0.1)) {
    plot <- plot + annotate("text",x=stats$BioFraction,y=max(stats$ypos),label=stats$symbol,size=7)
  }
  # Add Custom colors and modify legend title and labels.
  mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
  plot <- plot + scale_colour_manual(name="Replicate", values=colors, labels=mylabs) 
  #plot <- plot + scale_x_discrete(labels=stats$Cfg.Force)
  #plot <- plot + xlab("Force (xg)")
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  # Edit y axis.
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  # Make x and y axes consistent.
  plot <- plot + theme(axis.text.x = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
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
} #EOF


## Set-up the workspace -------------------------------------------------------

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
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(msstats_results)

# combine protein data and statistical results
prot_df <- left_join(msstats_prot,msstats_results,by=c("Protein","BioFraction"))

# Colors for the plot.
# Generated shade of black online: https://mycolor.space/
graphpad_purple <- c("R"=148,"G"=33,"B"=146)
colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
	    "#942192","#B847B4","#DC6AD7") # Swip Purples
	

## Generate the plots ---------------------------------------------------------

# sort proteins by module membership, drop M0
sorted_prots <- unlist(split(names(partition),partition)[-1])

# Loop to generate plots for all_proteins.
message("\nGenerating plots for ", 
	formatC(length(sorted_prots),big.mark=","), " proteins.")
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
	plots[[protein]] <- plot_protein(prot_df,gene_map,protein)
	setTxtProgressBar(pbar,value=match(protein,sorted_prots))
}
close(pbar)


# Generate a legend ------------------------------------------------------------

# Generate a plot with a legend.
plot <- plot_protein(prot_df,gene_map,swip, legend = TRUE)
plot_legend <- cowplot::get_legend(plot) 
myfile <- file.path(figsdir,"Protein_plots_legend.pdf")
ggsave(plot_legend,file=myfile,width=4.5,height=4.5)


## Annotate plots with module assignment ---------------------------------------

for (i in c(1:length(plots))) {
	protein <- names(plots)[i]
	plot <- plots[[protein]]
	module <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
		select(Abundance) %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plot <- plot + annotate(geom="label",x=7, y=ypos, label=module)
	plots[[protein]] <- plot
}


## Save the data --------------------------------------------------------------

# Save all proteins.
message("\nSaving all plots, this will take several minutes.")
myfile <- file.path(figsdir,"Protein_plots.pdf")
ggsavePDF(plots, myfile)
