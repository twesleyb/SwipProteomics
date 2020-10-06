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
# FIXME: whats required?
#devtools::load_all()

# other imports
suppressPackageStartupMessages({
	library(dplyr) 
	library(MSstatsTMT) # my fork
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
#myfile <- file.path(root,"rdata","msstats_prot.rda")
#load(myfile) # == msstats_prot
data(msstats_prot)

# munge - clarify covariate names
#conditition <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
#biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
mixture <- gsub("_1","",as.character(msstats_prot$Run))
#msstats_prot$Condition <- conditition
#msstats_prot$BioFraction <- biofraction
msstats_prot$Mixture <- mixture


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
#load(file.path(root,"rdata","msstats_contrasts.rda"))
data(msstats_contrasts.rda)
cm0 <- msstats_contrasts # intrafraction contrasts

knitr::kable(cm0[1,])


## fit model for swip --------------------------------------------------------
#fx <- formula("Abundance ~ 1 + (1|Mixture) + Condition * BioFraction")
#fx <- formula("Abundance ~ 1 + (1|Mixture) + Condition.BioFraction")
# NOTE: declaring the interaction of Condition.Biofraction like : doesnt work!
#fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
# NOTE: changed Run annotation to make more sense
# Avoid changing Condition as downstream functions seem to depened upon
# formatting of stuff

prot <- "Q3UMB9" # Swip
fx <- formula("Abundance ~ 1 + (1|Mixture) + Condition") # | condition is interaction(BioF,Cond)

fit <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(fit_list)

df <- bind_rows(results)

knitr::kable(df)

print(summary(fit))
