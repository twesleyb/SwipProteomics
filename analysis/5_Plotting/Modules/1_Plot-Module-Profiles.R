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
data(module_gof) # module gof stats!
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

# we will plot protein abundance adjusted for Batch effect (Mixture)
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

plot_profile <- function(module, msstats_prot, partition,
			 module_colors, module_gof, wt_color = "#47b2a4") {
	if (module %notin% module_gof$Module) {
	       warning(module," is not in 'module_gof'.")
	       return(NULL)
	}
  # Subset
  subdat <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein])) %>% 
	  filter(Module == module)
  # number of proteins in module
  nprots <- length(unique(subdat$Protein))
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
  # calculate coefficient of variation (CV == unitless error) and scale to max
  df <- df %>% mutate(CV = SD/mean_Abundance)
  df <- df %>% group_by(Protein) %>%
	mutate(scale_Abundance = mean_Abundance/max(mean_Abundance))
  # get module fitted data by fitting linear model to scaled Abundance
  # we do it this way because we want to plot the data w/o batch effect
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
  # get r2 for annot plot title
  # need to regen stats
  r2 <- module_gof %>% filter(Module == module) %>% select(R2.fixef) %>% as.numeric() 
  r2_anno <- paste("(",paste(paste(c("R2.Fixef = "),round(r2,3)),collapse=" | "),")")
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
  plot <- plot + ggtitle(paste0(module," (n = ",nprots,")\n",r2_anno))
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

#
#data(washc_prots)
#
#part <- setNames(rep(partition[swip],length(washc_prots)), nm = washc_prots)
#plot <- plot_profile("M17",msstats_prot,part,module_colors,module_gof)
#myfile <- file.path(figsdir, "washc_prots.pdf")
#ggsave(myfile, plot)


## generate plots -------------------------------------------------------------

# loop to generate plots for all modules
modules <- split(partition,partition)[-1]
names(modules) <- paste0("M",names(modules))

stopifnot(length(module_colors) >= length(modules)) # insufficient colors

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to generate plots
message("\nGenerating profile plots of ", length(modules), " modules.")

plot_list <- foreach(module = names(modules)) %dopar% {
	plot_profile(module, msstats_prot, partition, module_colors, module_gof)
}
names(plot_list) <- names(modules)

# drop null
bad_modules <- names(which(sapply(plot_list,is.null)))
idx <- names(plot_list) %notin% bad_modules

# save plots as a single pdf
message("\nSaving plots as a single pdf.")
myfile <- file.path(figsdir,"Module_Profiles.pdf")
ggsavePDF(plot_list[idx],myfile)
