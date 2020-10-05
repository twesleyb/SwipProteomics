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

# load projects functions in root/R
devtools::load_all()
# FIXME: whats required?

# other imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
	library(MSstatsTMT) # my fork
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
load(file.path(root,"rdata","msstats_contrasts.rda"))
cm0 <- msstats_contrasts # intrafraction contrasts


## fit model for swip --------------------------------------------------------

prot <- "Q3UMB9" # Swip
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
fit <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(fit_list)

df <- bind_rows(results)

knitr::kable(df)

print(summary(fit))
