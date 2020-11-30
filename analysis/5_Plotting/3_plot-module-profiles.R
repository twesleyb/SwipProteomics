#!/usr/bin/env Rscript

# title: SwipProteomics
# description: generate module-level profile plot
# author: Tyler W Bradshaw

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"

input_colors = "ne_surprise_colors"
input_part = "ne_surprise_partition"


## ---- Prepare the R environment

renv::load(root,quiet=TRUE)

devtools::load_all(root, quiet=TRUE)

# load the data
data(msstats_prot)
data(list=input_part) # partition
data(list=input_colors) # module_colors

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
	library(doParallel)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set plotting theme and font
ggtheme()
set_font("Arial", font_path=fontdir)


## ---- function 

module = "M1"
prots = names(partition[partition==1])

plotModule <- function(module, prots, msstats_prot, module_colors) {

  # color for Control condition
  wt_color = "#47b2a4"

  # Subset
  subdat <- msstats_prot %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # prepare the data
  df <- subdat %>% group_by(Protein) %>%
	  mutate(Intensity = 2^Abundance) %>% 
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Intensity = log2(median(rel_Intensity)),
	          SD = sd(log2(rel_Intensity)),
	          N = length(rel_Intensity),
	          .groups="drop") %>%
	  mutate(CV = SD/med_Intensity)

  df$scale_Intensity <- scale01(df$med_Intensity)

  # get module fitted data by fitting linear model to scaled Intensity
  # NOTE: no Mixture
  fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)
  fm <- lmerTest::lmer(fx, df)

  # collect coefficients
  fit_df <- data.table("coef" = names(lme4::fixef(fm)),
		       "fit_y" = lme4::fixef(fm)) %>%
    mutate(Genotype = gsub("Genotype","",sapply(strsplit(coef,"\\:"),"[",1))) %>%
    mutate(BioFraction=gsub("BioFraction","",sapply(strsplit(coef,"\\:"),"[",2)))

  # combine module data and fitted values
  df <- left_join(df, fit_df,by=c("Genotype","BioFraction"))

  # again, insure factor order is correct
  df$Genotype <- factor(df$Genotype,levels=c("Control","Mutant"))
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  
  # get marginal r2 for annot plot title
  r2 <- r.squaredGLMM.merMod(fm)
  r2_anno <- paste("(",paste(paste(c("R2.Fixef = "),
			     round(r2[,"R2m"],3)),collapse=" | "),")")

  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Intensity)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + aes(ymin=scale_Intensity - CV)
  plot <- plot + aes(ymax=scale_Intensity + CV)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(paste0(module," (n = ",nprots,")\n",r2_anno))
  plot <- plot + ylab("Scaled Protein Intensity")
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
  mut_color <- module_colors[[module]]

  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))

  return(plot)
} #EOF


## ---- generate plots 

# all modules
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

stopifnot(all(names(modules) %in% names(module_colors)))

message("\nGenerating profile plots of ", length(modules), " modules.")


# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to generate plots
plot_list <- foreach(module = names(modules)) %dopar% {
	plotModule(module, prots=modules[[module]], 
			     msstats_prot, module_colors)
} #EOL
names(plot_list) <- names(modules)

# drop null
bad_modules <- names(which(sapply(plot_list,is.null)))
idx <- names(plot_list) %notin% bad_modules


## ---- save plots as a single pdf

message("\nSaving plots as a single pdf.")

myfile <- file.path(figsdir,"module_profiles.pdf")
ggsavePDF(plot_list[idx],myfile)
