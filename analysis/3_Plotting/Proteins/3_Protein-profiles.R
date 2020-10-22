#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:

# Input data in root/data/
root = "~/projects/SwipProteomics"

## Prepare environment --------------------------------------------------------

renv::load(root,quiet=TRUE)

suppressWarnings({ devtools::load_all() })

# load the data
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# other imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs","Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

ggtheme(); set_font("Arial",font_path=fontdir)


## plot protein summary --------------------------------------------------------


plot_protein_summary <- function(protein, wt_color = "#47b2a4",
				 mut_color = col2hex(c("R"=148,"G"=33,"B"=146))) {
  # assumes msstats_prot and msstats_results are available in the current env
  # get protein's gene symbol
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  # annotate msstats_results with BioFraction
  biofraction <- sapply(strsplit(msstats_results$Label,"\\."),"[",3)
  prot_df <- msstats_prot
  results_df <- msstats_results
  results_df$BioFraction <- factor(biofraction,levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # prepare the data -- why not do it in one long pipe?
  df <- prot_df %>% left_join(results_df,by=c("Protein","BioFraction")) %>%
	  filter(Protein %in% names(partition)) %>%  
	  mutate(Module = paste0("M",partition[Protein])) %>%
	  filter(Protein == protein) %>% 
	  select(Protein, Mixture, Channel, Condition, BioFraction, Abundance, SE, Module) %>% 
	  unique() %>% 
	  group_by(Condition,Protein) %>% 
	  summarize(mean_Abundance = mean(Abundance),
		    SD = sd(Abundance),
		    SE = unique(SE),.groups="drop") %>% 
	  mutate(CV = SD/mean_Abundance) %>% 
	  mutate(norm_Abundance = mean_Abundance/max(mean_Abundance)) %>%
	  mutate(BioFraction = sapply(strsplit(as.character(Condition),"\\."),"[",2)) %>%
	  mutate(Genotype = sapply(strsplit(as.character(Condition),"\\."),"[",1)) %>%
	  mutate(BioFraction = factor(BioFraction, levels=c("F4","F5","F6","F7","F8","F9","F10")))
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = Genotype)
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Normalized Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.x = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.x = element_text(family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.y = element_text(family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  return(plot)
} #EOF


## generate plots -------------------------------------------------------------

# protein names sorted by module membership
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

# loop for all proteins
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
	module <- paste0("M",partition[protein])
	plot <- plot_protein_summary(protein,mut_color=module_colors[module])
	plot_label <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
		select(norm_Abundance) %>% range()
        ypos <- yrange[1] - 0.1* diff(yrange)
        plot <- plot + annotate(geom="label",x=7, y=ypos, label=plot_label)
        plots[[protein]] <- plot
        setTxtProgressBar(pbar,value=match(protein,sorted_prots))
} # EOL for proteins
close(pbar)

# save
myfile = file.path(figsdir,"Protein_profiles.pdf")
ggsavePDF(plots, file=myfile)
