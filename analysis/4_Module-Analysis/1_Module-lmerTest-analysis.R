#!/usr/bin/env Rscript

# title: SwipProteomics
# description: module-level analysis with mixed models
# author: Tyler W Bradshaw

## ---- Inputs

# Input data in root/data/
root = "~/projects/SwipProteomics"

input_part = "ne_surprise_partition"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load the data
data(gene_map) # gene_map
data(msstats_prot) # msstats_prot
data(list=input_part) # partition

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(doParallel)
})


## ---- Function


fitModule <- function(prots, msstats_prot) {
  # fitting mixed-model to log2 relative Intensity
  fx <- "log2(rel_Intensity) ~ 1 + Condition + (1|Protein) + (1|Mixture)"
  # build list of input args for lmerTest
  lmer_args <- list()
  lmer_args[["formula"]] <- fx
  # prepare the data -- scaled (relative) Intensity
  lmer_args[["data"]] <- msstats_prot %>% 
	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein) %>% 
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  subset(Protein %in% prots)
  fm <- do.call(lmerTest::lmer, lmer_args)
  # assess overall contrast and collect results
  LT <- getContrast(fm,"Mutant","Control")
  result <- lmerTestContrast(fm,LT) %>% 
	  mutate(Contrast='Mutant-Control') %>% unique() %>% 
	  mutate('nProts'=length(prots)) 
  return(result)
} #EOF


## ---- main

# loop to fit module-level models and assess overall contrast
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

message("k Modules: ", length(modules))


## ---- loop to fit module-level models and assess contrast

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

results_list <- foreach(module = names(modules)) %dopar% {
  fitModule(modules[[module]], msstats_prot)
} # EOL
names(results_list) <- names(modules)

## collect results
results_df = bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>% 
	arrange(Pvalue)

# nSig modules
m <- results_df %>% filter(Padjust < 0.05) %>% select(Module) %>% unlist(use.names=FALSE)

message("n Sig Modules: ", length(m), " (Padjust < 0.05)")


## ---- save results

# save sig modules
sig_modules <- m
myfile <- file.path(root,"data","sig_modules.rda")
save(sig_modules, file=myfile,version=2)

# write results to excel 
res_file <- file.path(root,"tables","SWIP_Module_Results.xlsx")
write_excel(list("Module Results" = results_df), res_file)

# save results as rda
module_results <- results_df
myfile <- file.path(root,"data","module_results.rda")
save(module_results, file=myfile, version=2)
