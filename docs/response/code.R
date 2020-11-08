#!/usr/bin/env Rscript

## Illustrate the Protein-level analysis

root <- "~/projects/SwipProteomics"
renv::load(root)

#<<fit the model, eval=TRUE, echo=TRUE>>=

library(dplyr)
library(data.table)


## load SwipProteomics data
# devtools::install_github("twesleyb/SwipProteomics")
library(SwipProteomics)

data(swip)
data(msstats_prot)

## formula to be fit:
fx0 <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)")

# fit the model for WASHC4
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

# create a contrast
coeff <- lme4::fixef(fm0)
contrast0 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast0["GenotypeMutant:BioFractionF7"] <- +1 # positive coeff
contrast0["GenotypeControl:BioFractionF7"] <- -1 # negative coeff

# evaluate contrast
lmerTestContrast(fm0,contrast0) %>% knitr::kable()

#@

#<<fit the model, eval=TRUE, echo=TRUE>>=

## assess 'Mutant-Control' comparison

## create a contrast
contrast1 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast1[grepl("Mutant",names(contrast1))] <- +1/7
contrast1[grepl("Control",names(contrast1))] <- -1/7

lmerTestContrast(fm0, contrast1) %>% mutate(Contrast='Mutant-Control') %>% 
	unique() %>% knitr::kable()

#@

#<<fit the model, eval=TRUE, echo=TRUE>>=

## create a contrast
contrast1 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast1[grepl("Mutant",names(contrast1))] <- +1/7
contrast1[grepl("Control",names(contrast1))] <- -1/7

lmerTestContrast(fm0, contrast1) %>% mutate(Contrast='Mutant-Control') %>% 
	unique() %>% knitr::kable()

#@

#<<fit the model, eval=TRUE, echo=TRUE>>=

## illustrate module-level analysis

data(washc_prots) # "Q8C2E7" "Q6PGL7" "Q3UMB9" "Q9CR27" "Q8VDD8"

# the module-level formula:
fx1 <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)")

# fit the model
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein %in% washc_prots))

# assess contrast
lmerTestContrast(fm1, contrast1) %>% mutate(Contrast='Mutant-Control') %>% 
	unique() %>% knitr::kable()

#@


#<<fit the model, eval=TRUE, echo=TRUE>>=

## assess gof with Nakagawa coefficient of determination
r.squaredGLMM.merMod(fm1) %>% knitr::kable()


## calculate the variance components with variancePartition

# variancePartition expects all factors to be mixed effects
form <- formula(Abundance ~ (1|Genotype) + (1|BioFraction) + 
		(1|Mixture) + (1|Protein))

fit <- lme4::lmer(form, msstats_prot %>% filter(Protein %in% washc_prots))

variancePartition::calcVarPart(fit)

#@
