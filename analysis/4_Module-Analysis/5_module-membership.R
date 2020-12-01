#!/usr/bin/env Rscript

# title: SwipProteomics
# author: Tyler W Bradshaw
# description:

## ---- Inputs

# Input data in root/data/
root = "~/projects/SwipProteomics"
input_part = "ne_surprise_partition"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load the data
data(msstats_prot) # msstats_prot
data(list=input_part) # partition

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(doParallel)
})


## ---- function 

moduleMembership <- function(prots, msstats_prot) {
  # lmer based module membership
  dm <- msstats_prot %>%
	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein) %>%
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  group_by(Protein, Condition) %>%
	  summarize(scale_Intensity = log2(median(rel_Intensity)),
	  	   .groups="drop") %>% 
	  reshape2::dcast(Protein ~ Condition, value.var = "scale_Intensity") %>%
	  as.data.table() %>%
	  as.matrix(rownames="Protein")
  # fitting mixed-model to log2 relative Intensity
  fx <- scale_Intensity ~ 0 + Condition + (1|Protein)
  # build list of input args for lmerTest
  lmer_args <- list()
  lmer_args[["formula"]] <- fx
  lmer_args[["data"]] <- msstats_prot %>% 
	  subset(Protein %in% prots) %>%
	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein) %>%
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  group_by(Protein, Condition) %>%
	  summarize(scale_Intensity = log2(median(rel_Intensity)),
	  	  .groups="drop") 
  # fit the model with some lmer control
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- do.call(lmerTest::lmer, lmer_args)
  # get fixef
  y <- lme4::fixef(fm)
  names(y) <- gsub("Condition","",names(y))
  # examine cor between proteins and module fit y
  mm <- cor(t(dm[prots,]),y, method="pearson")
  kme <- setNames(as.numeric(mm),rownames(mm))
  return(kme)
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
  moduleMembership(modules[[module]], msstats_prot)
} # EOL
names(results_list) <- names(modules)

## collect results
module_membership = results_list


## ---- save results

myfile <- file.path(root, "data", "module_membership.rda")
save(module_membership,file=myfile,version=2)
