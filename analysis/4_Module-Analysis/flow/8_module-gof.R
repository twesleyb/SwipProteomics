#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(swip)
data(gene_map)
data(sigprots)
data(partition)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(variancePartition)
})


## functions ------------------------------------------------------------------

moduleGOF <- function(module,partition,msstats_prot){
  # a function the evaluates gof of modules with variancePartition
  # NOTE: variancePartition expects all factors to be modeled as mixed effects
  prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein]))
  form1 <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Protein)
  # fit the model for variancePartition
  fm1 <- lmer(form1, prot_df %>% filter(Module==module))
  vpart <- variancePartition::calcVarPart(fm1)
  # fit the model for calculating Nakagawa R2
  form2 <- Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)
  fm2 <- lmer(form2, prot_df %>% filter(Module==module))
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm2)),
		 nm=c("R2.fixef","R2.total"))
  rho <- c(vpart,r2)
  # fit 1 is the model for variancePartition; fit2 is the fit used for stats
  rho[["vp.isSingular"]] <- lme4::isSingular(fm1)
  rho[["lmer.isSingular"]] <- lme4::isSingular(fm2)
  return(rho) # gof stats
} #EOF


## fit all module models -------------------------------------------------------

# analysis of intra-module variance
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# loop to evaluate gof
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

# drop null results
idx <- sapply(results_list,is.null)
message("There were problems fitting ", sum(idx), " models.")

# collect results
module_sizes <- sapply(modules,length)
df <- as.data.table(do.call(rbind,results_list[!idx]),keep.rownames="Module")
df <- df %>% arrange(desc(Genotype))
df <- df %>% mutate(Size = module_sizes[Module])
df$"vp.isSingular" <- NULL
df$"lmer.isSingular" <- NULL
module_gof <- df %>% 
	select(Module, Size, BioFraction, Genotype,
	       Mixture, Protein, Residuals, R2.fixef, R2.total)

# work
module_gof %>% mutate(Quality = R2.total * sum(BioFraction, Genotype)/sum(Protein,Residuals))

q <- sum(module_gof$Quality)/length(modules)
message("Partition Quality: ", round(q,5))

module_gof %>% knitr::kable()


## save results -----------------------------------------------------------------

# save as rda
myfile <- file.path(root,"data","module_gof.rda")
save(module_gof, file=myfile, version=2)

# save the data
fwrite(module_gof,file.path(root,"rdata","module_gof.csv"))
