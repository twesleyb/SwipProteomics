#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:

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


plot_protein <- function(prot_df, gene_map, protein, 
			 sigprots, protein_gof, legend=FALSE) {
  # a function that generates the plot
  #wt_color = "#47b2a4"
  #mut_color <- col2hex(c("R"=148,"G"=33,"B"=146))
  title_colors <- c("red"=TRUE,"black"=FALSE)
  colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
  	    "#942192","#B847B4","#DC6AD7") # Swip Purples
  r2 <- protein_gof %>% filter(Protein == protein) %>% 
	select(R2.total) %>% as.numeric()
  title_anno <- paste0("(R2 = ",round(r2,3),")")
  # function to generate a proteins plot
  # Subset the data.
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  title_color <- names(which(title_colors == (protein %in% sigprots)))
  df <- subset(prot_df,prot_df$Protein == protein)
  # Insure Fraction is a factor, and levels are in correct order.
  df$BioFraction <- factor(df$BioFraction,
  			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Collect FDR stats.
  stats <- df %>% group_by(Genotype,BioFraction) %>% 
  		summarize(`Max Abundance` = max(norm_Abundance),
  			  FDR = unique(FDR),.groups="drop")
  stats$ypos <- 1.02 * max(stats$`Max Abundance`)
  stats <- stats %>% filter(Genotype == "Control")
  stats$symbol <- ""
  stats$symbol[stats$FDR<0.1] <- "."
  stats$symbol[stats$FDR<0.05] <- "*"
  stats$symbol[stats$FDR<0.005] <- "**"
  stats$symbol[stats$FDR<0.0005] <- "***"
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction, y = norm_Abundance)
  plot <- plot + aes(group = interaction(Mixture,Genotype))
  plot <- plot + aes(colour = interaction(Mixture,Genotype))
  plot <- plot + aes(shape=Mixture)
  #plot <- plot + aes(fill=Genotype)
  plot <- plot + aes(group = interaction(Mixture,Genotype))
  plot <- plot + aes(colour = interaction(Mixture,Genotype))
  plot <- plot + geom_point(size=2)
  plot <- plot + geom_line()
  plot <- plot + ggtitle(paste(gene,"|",protein,title_anno))
  plot <- plot + theme(plot.title=element_text(color=title_color))
  # Annotate with significance stars.
  if (any(stats$FDR<0.1)) {
    plot <- plot + annotate("text", 
			    x=stats$BioFraction, 
			    y=max(stats$ypos), 
			    label=stats$symbol,size=7)
  }
  # Add Custom colors and modify legend title and labels.
  mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
  plot <- plot + scale_colour_manual(name="Subject", values=colors,labels=mylabs) 
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

# load renv
root <- getrd()
renv::load(root,quiet=TRUE)

# Load additional functions in root/R/
devtools::load_all()

# load the data
data(swip)
data(gene_map)
data(partition)
data(protein_gof)
data(msstats_prot)
data(module_colors)
data(msstats_results)

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})


# Project directories:
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# Set theme for the plots:
ggtheme(); set_font("Arial",font_path=fontdir)


# combine protein data and statistical results
# we will use stats to annotate plots with stars
biofraction <- sapply(strsplit(msstats_results$Contrast,"\\."),"[",3)
msstats_results$BioFraction <- biofraction
shared_cols <- intersect(colnames(msstats_prot),colnames(msstats_results))
prot_df <- left_join(msstats_prot,msstats_results,by=shared_cols)

# annotate with module membership
prot_dt <- prot_df %>% filter(Protein %in% names(partition))
prot_df$Module <- paste0("M",partition[prot_df$Protein])

# collect vector of sigprots. annotate title in red if protein has overall 
# sig change in 'Mutant-Control' comparison
sigprots <- msstats_results %>% filter(Contrast == 'Mutant-Control') %>% 
	filter(FDR<0.05) %>% select(Protein) %>% 
	unlist() %>% as.character() %>% unique()


## Generate the plots ---------------------------------------------------------

# sort proteins by module membership, drop M0
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

# Loop to generate plots for all_proteins.
message("\nGenerating plots for ", 
	formatC(length(sorted_prots),big.mark=","), " proteins.")

plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)

for (protein in sorted_prots) {
	# generate a proteins plot
	plot <- plot_protein(prot_df,gene_map,protein,sigprots,protein_gof)
	# annotate with module assignment
	plot_label <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
		select(norm_Abundance) %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plots[[protein]] <- plot + annotate(geom="label",x=7, y=ypos, 
					    label=plot_label)
	setTxtProgressBar(pbar,value=match(protein,sorted_prots))
}
close(pbar)


# Generate a plot with a legend.
plot <- plot_protein(prot_df,gene_map, swip, sigprots, protein_gof, legend = TRUE)
plot_legend <- cowplot::get_legend(plot) 


## save ------------------------------------------------------------

# legend
myfile <- file.path(figsdir,"Protein_plots_legend.pdf")
ggsave(plot_legend,file=myfile,width=4.5,height=4.5)

# plot list
message("\nSaving ",length(plots), " plots as a single PDF.")
myfile <- file.path(figsdir,"Protein_plots.pdf")
ggsavePDF(plots,myfile)
