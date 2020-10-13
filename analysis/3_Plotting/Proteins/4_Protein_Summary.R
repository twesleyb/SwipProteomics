#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
wt_color = "#47b2a4"
mut_color <- c("R"=148,"G"=33,"B"=146)

# Input data in root/data/
root = "~/projects/SwipProteomics"


## Output ---------------------------------------------------------------------


## Functions -------------------------------------------------------------------


## Prepare environment --------------------------------------------------------

renv::load(root,quiet=TRUE)

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

suppressWarnings({ devtools::load_all() })

fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

ggtheme(); set_font("Arial",font_path=fontdir)

# load the data
data(gene_map)
data(partition)
data(sig_modules)
data(msstats_prot)
data(module_colors)
data(msstats_results)

## plot protein summary --------------------------------------------------------

plot_protein_summary <- function(protein, wt_color, mut_color) {
  # a function to generate a proteins summary plot
  require(dplyr,quietly=TRUE)
  require(ggplot2,quietly=TRUE)
  msstats_df <- left_join(msstats_prot,msstats_results,by=c("Protein","BioFraction"))
  msstats_df$Module <- paste0("M",partition[msstats_df$Protein])
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  df <- subset(msstats_df,msstats_df$Protein == protein) %>% 
	select(Mixture,Channel,Condition,Protein,Abundance,SE,Module) %>% 
	unique() %>% 
	group_by(Condition,Protein) %>% 
	summarize(mean_Abundance = mean(Abundance),
		  SD = sd(Abundance),
		  SE = unique(SE),.groups="drop")
  df <- df %>% mutate(CV = SD/mean_Abundance)
  df <- df %>% mutate(norm_Abundance = mean_Abundance/max(mean_Abundance))
  df$BioFraction <- sapply(strsplit(as.character(df$Condition),"\\."),"[",2)
  df$Genotype <- sapply(strsplit(as.character(df$Condition),"\\."),"[",1)
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = Genotype, colour = Genotype, 
		     shape=Genotype, fill=Genotype,shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Normalized Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black", size=11, 
						  angle = 0, hjust = 1, 
						  family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",size=11, 
						  angle = 0, hjust = 1, 
						  family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  return(plot)
} #EOF


## ----------------------------------------------------------------------------

# Loop to do work.
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

for (module in names(modules)) {
  plots <- list()
  message("\nWorking on: ", module)
  proteins <- modules[[module]]
  pbar <- txtProgressBar(max=length(proteins),style=3)
  for (protein in proteins) {
    # FIXME: colors do not seem correct
    plot <- plot_protein_summary(protein,wt_color,col2hex(mut_color))
    plot_label <- paste("Module:", partition[protein])
    yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
  	    select(norm_Abundance) %>% range()
    ypos <- yrange[1] - 0.1* diff(yrange)
    plot <- plot + annotate(geom="label",x=7, y=ypos, label=plot_label)
    plots[[protein]] <- plot
    setTxtProgressBar(pbar,value=match(protein,proteins))
  } # EOL for proteins
  close(pbar)
  # save
  myfile = file.path(figsdir,"Modules",paste0(module,"_Protein_summary.pdf"))
  ggsavePDF(plots,file=myfile)
  message("\nSaved: ", module)
} # EOL for modules

#fx = "Abundance ~ 0 + (1|BioFraction) + Genotype"
#fm = lmerTest::lmer(fx, msstats_prot %>% filter(Protein == swip))
