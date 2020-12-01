#!/usr/bin/env Rscript

# title: SwipProteomics
# description: generate module-level profile plot
# author: Tyler W Bradshaw

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

input_colors = "ne_surprise_colors"
input_part = "ne_surprise_partition"


## ---- Prepare the R environment

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load the data
data(module_gof)
data(sig_modules)
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


plotModule <- function(module, prots, msstats_prot, 
		       module_colors, module_gof, title_color="black") {

  # color for Control condition
  wt_color = "#47b2a4"

  # plot title annot
  #r2 <- module_gof %>% filter(Module == module) %>% 
  #	  select(R2.fixef, R2.total) %>% as.numeric()
  #title_anno <- paste(paste0(c("(R2_fixef = ","(R2_fixef = "), 
  #			     round(r2,3),")"),collapse=" ")

  # Subset
  subdat <- msstats_prot %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # prepare the data
  df <- subdat %>% 
	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein) %>%
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Intensity = median(rel_Intensity),
	          .groups="drop") %>%
	  mutate(scale_Intensity = scale01(log2(med_Intensity)))

  # get module fitted data by fitting linear model to scaled Intensity
  # explicitly estimate all coeff by setting intercept to 0
  fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, df, control=lmer_control)

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

  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Intensity)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
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
  plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
			   linetype="dashed",alpha=1,size=0.75)
  mut_color <- module_colors[[module]]
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  #plot <- plot + ggtitle(paste0(module," (n = ",nprots,")\n",title_anno))
  plot <- plot + ggtitle(paste0(module," (n = ",nprots,")"))
  plot <- plot + theme(plot.title = element_text(color=title_color))

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
	plotModule(module, prots=modules[[module]], msstats_prot, 
		   module_colors, module_gof)
} #EOL
names(plot_list) <- names(modules)


## ---- save plots as a single pdf

message("\nSaving plots as a single pdf.")

myfile <- file.path(figsdir,"module_profiles.pdf")
ggsavePDF(plot_list,myfile)
