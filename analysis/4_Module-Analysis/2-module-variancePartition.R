#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

input_part <- "ne_surprise_surprise_partition"

save_results = TRUE

## ---- prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)


## ---- load the data
data(swip)
data(gene_map)
data(msstats_prot)
data(list=input_part)


## ---- requirements
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(variancePartition)
})


## ---- function

#module = "M1"
#table(partition)
#head(msstats_prot)

moduleGOF <- function(module, partition, msstats_prot){
  # a function the evaluates gof of a module with variancePartition
  # NOTE: variancePartition expects all factors to be modeled as mixed effects

  ## subset the data, annotate with module membership, calc scale_Intensity
  prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein])) %>% 
	  group_by(Protein) %>% 
	  mutate(Intensity = 2^Abundance) %>% 
	  mutate(scale_Intensity = Intensity/sum(Intensity))


  ## the formula for variancePartition -- all effects modeled as random effects
  form <- log2(scale_Intensity) ~ (1|Mixture) + (1|Genotype) + 
	  (1|BioFraction) + (1|Protein)

  # fit the model with some lmerControl
  lmer_control <- lme4::lmerControl(check.conv.singular = "ignore")
  vpart_fm <- lme4::lmer(form, data = prot_df %>% filter(Module==module), 
			 control=lmer_control)

  # do the variancePartition bit == vp=getVariance(fm); pve=vp/sum(vp)
  vpart <- variancePartition::calcVarPart(vpart_fm)


  ## fit the model for calculating Nakagawa R2
  fx <- log2(scale_Intensity) ~ 1 + Condition + (1|Protein)
  fm <- lme4::lmer(fx, prot_df %>% filter(Module==module), 
		    control = lmer_control)


  ## r.squaredGLMM.merMod
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),
		 nm=c("R2.fixef","R2.total"))

  # combine gof stats
  rho <- c(vpart, r2)

  rho[["isSingular"]] <- lme4::isSingular(fm)

  return(rho) # gof stats
} #EOF


## ---- fit all module-level models 

# analysis of intra-module variance
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# loop with pbar to evaluate gof for every module
message("\nEvaluating goodness-of-fit of modules.")
results_list <- list()
pbar <- txtProgressBar(max=length(modules),style=3)
for (module in names(modules)) {
	gof <- tryCatch(expr = { moduleGOF(module, partition, msstats_prot) },
			error = function(e) {}, # return null if error or 
			warning = function(w) {}) # warning
	results_list[[module]] <- gof
	setTxtProgressBar(pbar,value=match(module,names(modules)))
}
close(pbar)

# NOTE: the numerous messages about boundary/singular fits comes from the 
# from the variancePartition side of things. 
# From inspection this arises because Mixture often doesn't contribute much 
# to the module-level variance. This is not problematic.
# These warnings are suppressed with lmerControl.

# check for NULL results
idx <- sapply(results_list,is.null)
stopifnot(!any(idx))

# collect results from loop
module_sizes <- sapply(modules,length)
df <- as.data.table(do.call(rbind,results_list[!idx]),keep.rownames="Module")
df <- df %>% arrange(desc(Genotype))
df <- df %>% mutate(Size = module_sizes[Module])


## ---- save results 

if (save_results) {

  # save as rda
  module_gof <- df
  myfile <- file.path(root,"data","module_gof.rda")
  save(module_gof, file=myfile, version=2)
  
  # save the data
  results_list <- list("Module GOF"=module_gof)
  write_excel(results_list,file.path(root,"tables","S5_SWIP-TMT_Module_GOF.xlsx"))

}
