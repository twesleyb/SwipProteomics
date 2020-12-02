#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: 

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)

# load the data
data(gene_map)
data(select_paths)
data(msstats_prot)

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


## ---- module-level statistical analysis

# loop to do module-level lmer statistical analysis
results_list <- list()

pbar <- txtProgressBar(max=length(select_paths),style=3)
for (path in names(select_paths)){
  # build input args list for lmerTest
  lmer_args <- list()
  prots <- select_paths[[path]]
  # formula using scaled (relative) intensity
  # NOTE: intercept matters here!
  # NOTE: since this is a protein complex, maybe we can use Abundance and not
  # Intensity?
  fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)
  lmer_args[["formula"]] <- fx
  # subset the data
  lmer_args[["data"]] <- msstats_prot %>% 
	  subset(Protein %in% prots) %>% 
	  mutate(Intensity = 2^Abundance) %>%
	  group_by(Protein) %>% 
	  mutate(rel_Intensity=Intensity/sum(Intensity))
  # fit the model with lmer control
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore",
					      check.conv.grad="ignore")
  lmer_args[["control"]] <- lmer_control
  fm <- do.call(lmerTest::lmer, lmer_args)
  ## assess contrast
  LT <- getContrast(fm, "Mutant","Control")
  res <- lmerTestContrast(fm, LT) %>% 
	  mutate(Contrast="Mutant-Control") %>% unique()
  ## calc R2
  r2 <- r.squaredGLMM.merMod(fm)
  res$R2.fixef <- r2[,"R2m"]
  res$R2.total <- r2[,"R2c"]
  # other annot
  res$nProts <- length(prots)
  res$Pathway <- path
  ## evaluate variance explained by major covariates
  # NOTE: we ignore warnings about convergence and singular models
  # In these cases, a simpler model is prefered
  fx <- Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)
  fm <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% prots), 
		       control = lmer_control)
  vp <- getVariance(fm)
  pve <- vp/sum(vp)
  pve <- pve[names(pve) != "Fixed"] # fixef == 0
  names(pve) <- paste0("PVE_",names(pve))
  # return statistical results and variance, partitioned
  results_list[[path]] <- cbind(res, as.data.table(t(pve)))
  setTxtProgressBar(pbar,value=match(path,names(select_paths)))
}
close(pbar)

# FIXME: where is problem coming from?

## ---- function to plot a complex's profile

plotProfile <- function(path, prots, msstats_prot){

  # colors for Control and Mutant condition
  wt_color = "#47b2a4"
  mut_color = "#b671af"

  # subset
  subdat <- msstats_prot %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
  		    	     levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # prepare the data for plotting
  df <- subdat %>% 
  	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein) %>%
  	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
  	  group_by(Protein, Genotype, BioFraction) %>%
  	  summarize(med_Intensity = median(rel_Intensity), .groups="drop") %>%
  	  mutate(scale_Intensity = scale01(log2(med_Intensity/sum(med_Intensity))))

  # get module fit by fitting mixed model to scale Intensity
  fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)

  # fit the model with some lmer control
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, df, control = lmer_control)

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
  plot <- plot + ggtitle(paste0(path, " (n = ",nprots,")"))
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
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))

  return(plot)
} #EOF


## ---- loop to generate plots

plot_list <- list()
for (path in names(select_paths)){
	prots <- select_paths[[path]]
	plot_list[[path]] <- plotProfile(path, prots, msstats_prot)
} #EOL


## ---- save results

#myfile <- file.path(root,"figs","Modules","split-complexes.pdf")
myfile = "select_complexes.pdf"
ggsavePDF(plot_list, myfile)

myfile <- "complexes.csv"
dplyr::bind_rows(results_list) %>% fwrite(myfile)
