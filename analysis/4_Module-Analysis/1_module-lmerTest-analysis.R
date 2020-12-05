#!/usr/bin/env Rscript

# title: SwipProteomics
# description: module-level analysis with mixed models
# author: Tyler W Bradshaw

## ---- Inputs

# Input data in root/data/
root = "~/projects/SwipProteomics"
input_part = "ne_surprise_partition"

save_results = TRUE


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

  # fit mixed-model to log2 relative (scaled to sum) Intensity
  fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

  # build list of input args for lmerTest
  lmer_args <- list()
  lmer_args[["formula"]] <- fx
  lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots) %>% 
	  group_by(Protein) %>% mutate(Intensity = 2^Abundance) %>%
	  mutate(rel_Intensity=Intensity/sum(Intensity))

  # fit the model with some lmer control
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
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

# loop to do module-level analysis
results_list <- foreach(module = names(modules)) %dopar% {
  fitModule(modules[[module]], msstats_prot)
} # EOL
names(results_list) <- names(modules)

## collect results
results_df = bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>% 
	arrange(Pvalue)


## ---- save results

if (save_results) {

  # write results to excel 
  
  # drop singular col 
  results_df$isSingular <- NULL
  
  # re-arrange column order
  results_df <- results_df %>% 
  	dplyr::select(Module, nProts, Contrast, log2FC, 
			percentControl, SE, Tstatistic, 
			Pvalue, FDR, Padjust, DF, S2)
  
  # annotate candidate sig modules
  results_df <- results_df %>% 
	  mutate(candidate = Padjust < 0.05 & (percentControl > 1.05 | percentControl < 0.95))
  results_df <- results_df %>% arrange(desc(candidate))
  
  # summary
  message("n Sig modules: ", sum(results_df$candidate))

  # list of results
  results_list <- list()

  # data.frame describing network partition
  idx <- match(names(partition),gene_map$uniprot)
  df <-  data.table(UniProt = names(partition), 
  		 Entrez = gene_map$entrez[idx],
  		 Symbol = gene_map$symbol[idx],
  		 Membership = partition)

  results_list[["Partition"]] <- df %>% arrange(Membership)
  results_list[["Module Results"]] <- results_df 
  
  # save in root/tables
  namen <- gsub("partition","S4_SWIP-TMT_Module_Results.xlsx",input_part)
  myfile <- file.path(root,"tables",namen)
  write_excel(results_list, myfile)
  
  # save results as rda in root/data
  module_results <- results_df
  namen <- gsub("partition","module_results.rda",input_part) # e.g. ne_surprise2_module_results.rda
  myfile <- file.path(root,"data", namen)
  save(module_results, file=myfile, version=2)
 
}
