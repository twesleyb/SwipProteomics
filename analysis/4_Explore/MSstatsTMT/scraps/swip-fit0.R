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
knitr::kable(cm0[1,])


## fit model for swip --------------------------------------------------------

prot <- "Q3UMB9"

fx <- formula("Abundance ~ 1 + (1|Mixture) + Condition")

fit <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(fit_list)

df <- bind_rows(results)

knitr::kable(df)

print(summary(fit))

names(fit_list[[1]])

fit_list[[1]]$coeff

fit_list[[1]]$sigma

fit_list[[1]]$s2_df

fit_list[[1]]$s2

fit_list[[1]]$A
