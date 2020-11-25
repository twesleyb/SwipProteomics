#!/usr/bin/env Rscript

# title: SwipProteomics
# description:
# author: Tyler W Bradshaw

## ---- Inputs

# options:

# Input data in root/data/
root = "~/projects/SwipProteomics"


## ---- Prepare the R environment

renv::load(root,quiet=TRUE)

devtools::load_all(root,quiet=TRUE)

# load the data
data(partition)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(doParallel)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules")
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

## ---- Function


fit_profile <- function(prots, msstats_prot) {
  # Subset
  subdat <- msstats_prot %>% subset(Protein %in% prots)
  # number of proteins in module
  nprots <- length(unique(subdat$Protein))
  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # scale to max, take median of three replicates
  df <- subdat %>% group_by(Protein) %>%
	  mutate(scale_Abundance = Abundance/max(Abundance)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Abundance = median(scale_Abundance), 
	          SD = sd(scale_Abundance),
	          N = length(scale_Abundance),
	          .groups="drop")
  # calculate coefficient of variation (CV == unitless error) and scale to max
  df <- df %>% mutate(CV = SD/med_Abundance)
  # get module fitted data by fitting linear model to scaled Abundance
  fm <- lmerTest::lmer(med_Abundance ~ 0 + Genotype:BioFraction + (1|Protein), df)
  # assess overall contrast
  LT <- getContrast(fm,"Mutant","Control")
  res <- lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% unique()
  return(res)
} #EOF


## ---- main

# loop to generate plots for all modules
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to fit models, assess contrasts
results_list <- foreach(module = names(modules)) %dopar% {
  fit_profile(prots=modules[[module]], msstats_prot)
} # EOL
names(results_list) <- names(modules)

## collect results
results_df = bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>% 
	arrange(Pvalue)

## save
output_results <- file.path(root,"tables","SWIP_Module_Results.csv")
fwrite(results_df, output_results)
