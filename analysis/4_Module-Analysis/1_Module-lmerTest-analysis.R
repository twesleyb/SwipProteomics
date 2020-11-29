#!/usr/bin/env Rscript

# title: SwipProteomics
# description: module-level analysis with mixed models
# author: Tyler W Bradshaw

## ---- Inputs

# Input data in root/data/
root = "~/projects/SwipProteomics"

input_gene = "gene_map"
input_data = "msstats_prot"
input_part = "partition"

# * (gene_map)
# * (msstats_prot)
# * (partition)


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)

devtools::load_all(root, quiet=TRUE)

# load the data
data(list=input_data) # msstats_prot
data(list=input_part) # partition
data(list=input_gene) # gene_map

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(doParallel)
})


## ---- Function


fit_module <- function(prots, msstats_prot, gene_map) {

  # subset data for all prots in a module
  subdat <- msstats_prot %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  # always do this explicitly before plotting to save yourself trouble
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # for each protein: scale Abundance to max and take median of three* replicates
  #* NOTE: there may be less than three BioReplicates if a protein was identified in
  # less than 3 of the TMT mixtures. 
  df <- subdat %>% group_by(Protein) %>%
	  mutate(scale_Abundance = Abundance/max(Abundance)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Abundance = median(scale_Abundance), 
	          SD = sd(scale_Abundance),
	          N = length(scale_Abundance),
	          .groups="drop")

  # calculate coefficient of variation (CV == unitless error)
  df <- df %>% mutate(CV = SD/med_Abundance)

  # get module fitted data by fitting linear model to scaled Abundance
  fx <- "med_Abundance ~ 0 + Genotype:BioFraction + (1|Protein)"
  fm <- lmerTest::lmer(formula=fx, data=df)

  # assess overall contrast
  LT <- getContrast(fm,"Mutant","Control")
  result <- lmerTestContrast(fm,LT) %>% 
	  mutate(Contrast='Mutant-Control') %>% unique() 

  # annotate results with protein info
  idx <- match(prots, gene_map$uniprot)
  gene_prot <- paste(paste(gene_map$symbol[idx],prots,sep="|"),collapse="; ")
  result <- result %>% mutate("Proteins" = gene_prot) %>% mutate('nProts'=nprots) 

  return(result)
} #EOF


## ---- main

# loop to fit module-level models and assess overall contrast
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

message("N Modules: ", length(modules))

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

## loop to fit model and assess contrast
results_list <- foreach(module = names(modules)) %dopar% {
  fit_module(prots=modules[[module]], msstats_prot, gene_map)
} # EOL
names(results_list) <- names(modules)

## collect results
results_df = bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>% 
	arrange(Pvalue)


# nSig modules
m1 <- results_df %>% filter(FDR < 0.05) %>% select(Module) %>% unlist(use.names=FALSE)
m2 <- results_df %>% filter(Padjust < 0.05) %>% select(Module) %>% unlist(use.names=FALSE)
message("n Sig Modules: ", length(m1), " (FDR < 0.05)")
message("n Sig Modules: ", length(m2), " (Padjust < 0.05)")


## ---- save results

# save sig modules
sig_modules <- m2
myfile <- file.path(root,"data","sig_modules.rda")
save(sig_modules, file=myfile,version=2)

# write results to excel
output_results <- file.path(root,"tables","SWIP_Module_Results.xlsx")
write_excel(list("Module Results" = results_df), output_results)

# save results as rda
module_results <- results_df
myfile <- file.path(root,"data","module_results.rda")
save(module_results,file=myfile,version=2)
