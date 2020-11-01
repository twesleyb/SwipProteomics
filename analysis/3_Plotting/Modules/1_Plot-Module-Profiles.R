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

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
	library(doParallel)
})

stopifnot("norm_Abundance" %in% colnames(msstats_prot))

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules")
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set plotting theme and font
ggtheme(); set_font("Arial",font_path=fontdir)


## Function ------------------------------------------------------------------


# a function to generate protein profile plot:
plot_profile <- function(module, msstats_prot, partition,
			 module_colors, wt_color = "#47b2a4") {
  # Subset
  subdat <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein])) %>% 
	  filter(Module == module)
  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # calculate protein-wise mean of three replicates
  df <- subdat %>% group_by(Protein, Genotype, BioFraction) %>% 
        summarize(mean_Abundance = mean(norm_Abundance), # using norm_Abundance!
	          SD = sd(norm_Abundance),
	          N = length(norm_Abundance),
	          .groups="drop")
  # calculate coefficient of variation (CV == unitless error)
  df <- df %>% mutate(CV = SD/mean_Abundance)
  # scale protein profiles to maximum
  df <- df %>% group_by(Protein) %>%
	mutate(scale_Abundance = mean_Abundance/max(mean_Abundance))
  # get module fitted data by fitting linear model to scaled Abundance
  fm <- lm(scale_Abundance ~ 0 + Genotype:BioFraction, df)
  fit_df <- data.table("coef" = names(stats::coef(fm)),
		       "fit_y" = stats::coef(fm)) %>%
    mutate(Genotype = gsub("Genotype","",sapply(strsplit(coef,"\\:"),"[",1))) %>%
    mutate(BioFraction=gsub("BioFraction","",sapply(strsplit(coef,"\\:"),"[",2)))
  # combine module data and fitted values
  df <- left_join(df, fit_df,by=c("Genotype","BioFraction"))
  # again, insure factor order is correct
  df$Genotype <- factor(df$Genotype,levels=c("Control","Mutant"))
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Abundance)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + aes(ymin=scale_Abundance - CV)
  plot <- plot + aes(ymax=scale_Abundance + CV)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(module)
  plot <- plot + ylab("Scaled Abundance")
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
  # add fitted lines
  plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
			   linetype="dashed",alpha=1,size=0.75)
  # set colors
  plot <- plot + scale_colour_manual(values=c(wt_color,module_colors[[module]]))
  return(plot)
} #EOF


## generate plots -------------------------------------------------------------


# loop to generate plots for all modules
modules <- split(partition,partition)[-1]
names(modules) <- paste0("M",names(modules))

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to generate plots
plot_list <- foreach(module = names(modules)) %dopar% {
	plot_profile(module,msstats_prot,partition,module_colors)
}
names(plot_list) <- modules

# FIXME: annotate with lmer info and r2 and pve?

# save plots as a single pdf
myfile <- file.path(figsdir,"Module_Profiles.pdf")
ggsavePDF(plot_list,myfile)
