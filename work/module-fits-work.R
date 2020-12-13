#!/usr/bin/env Rscript

root <- getrd()
renv::load(root)
devtools::load_all(root)

library(dplyr)
library(data.table)

data(washc_prots)
data(msstats_prot)


fitModule <- function(fx) {
  args_list <- list()
  args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots) %>% 
	mutate(scale_Abundance = Abundance / max(Abundance))
  args_list[["formula"]] <- fx
  fm = lmerFit(args_list)
  LT = getContrast(fm,"Mutant","Control")
  lmerTestContrast(fm, LT) %>% mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()
}

fx = Abundance ~ 0 + Condition + Protein + (1|Mixture)
fitModule(fx)

fx1 = Abundance ~ 1 + Condition + Protein + (1|Mixture)
fitModule(fx1)

fx2 = Abundance ~ 1 + Condition + Protein + (1|Mixture)
fitModule(fx2)

fx3 = scale_Abundance ~ 1 + Condition + Protein + (1|Mixture)
fitModule(fx3)

fx4 = scale_Abundance ~ 1 + Condition + (1|Protein) + (1|Mixture)
fitModule(fx3)



# Condition Control is intercept, protein is fixed effect
fx = Abundance ~ Genotype + BioFraction + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots)
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
# create comparison to control
LT = fixef(fm)
LT[] <- 0 
LT["GenotypeMutant"] <- 1
calcSatterth(fm,LT)
lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()


##  how to define contrast when multiple fixef?
fx = Abundance ~ 1 + Condition + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots)
args_list[["formula"]] <- fx
fm = lmerFit(args_list)


## Protein as a fixed effect, intercept = 1

fx = Abundance ~ 1 + Condition + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots)
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()



## if scale, then logfc seems messed up
fx = scale_Abundance ~ 1 + Condition + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots) %>% group_by(Protein) %>% mutate(scale_Abundance =  Abundance/max(Abundance))
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()

# SCALE + Condition Control is intercept, protein is fixed effect
fx = scale_Abundance ~ Genotype + BioFraction + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots) %>% group_by(Protein) %>% mutate(scale_Abundance =  Abundance/max(Abundance))
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
# create comparison to control
LT = fixef(fm)
LT[] <- 0 
LT["GenotypeMutant"] <- 1
lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()

# FC estimate is still suppressed
# This probably means we are not defining contrast matrix correctly

fx = Abundance ~ 0 + Condition + Protein + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots) %>% group_by(Protein) %>% mutate(scale_Abundance =  Abundance/max(Abundance))
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
# create comparison to control
LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>% mutate(Contrast='Mutant-Control') %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()
r.squaredGLMM.merMod(fm) 


fx = Abundance ~ 0 + Condition + (1|Protein) + (1|Mixture)
args_list <- list()
args_list[["data"]] <- msstats_prot %>% subset(Protein %in% washc_prots) %>% group_by(Protein) %>% mutate(scale_Abundance =  Abundance/max(Abundance))
args_list[["formula"]] <- fx
fm = lmerFit(args_list)
# create comparison to control
LT = getContrast(fm,"ConditionMutant.F5","ConditionControl.F5")
lmerTestContrast(fm,LT) %>% mutate(Pvalue=formatC(Pvalue)) %>% unique() %>% knitr::kable()
