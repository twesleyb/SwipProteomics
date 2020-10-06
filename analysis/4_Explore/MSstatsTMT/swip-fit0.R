#!/usr/bin/env Rscript 

# title:
# author: twab
# description: Working through MSstatsTMT protein-level models and 
#   statistical comparisons


## prepare environment --------------------------------------------------------

# project root dir:
root <- "~/projects/SwipProteomics"

# load renv
renv::load(root)

# other imports
suppressPackageStartupMessages({
	library(dplyr) 
	library(MSstatsTMT) # my fork
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
devtools::load_all() 
data(msstats_prot)

# munge - clarify covariate names
mixture <- gsub("_1","",as.character(msstats_prot$Run))
msstats_prot$Mixture <- mixture


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
data(msstats_contrasts)
cm0 <- msstats_contrasts 

# all intrafraction contrasts in the format MSstatsTMT expects
# example, a single row/contrast:
#knitr::kable(cm0[1,])


## fit model for swip --------------------------------------------------------

prot <- "Q3UMB9"

fx0 <- formula("Abundance ~ 1 + (1|Mixture) + Condition")
message("#### fx0:",fx0)

fit <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx0, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)
# Error in as.formula(formula) : object 'fx' not found 

results <- MSstatsTMT::getResults(fit_list)

df <- bind_rows(results)

print(summary(fit))
