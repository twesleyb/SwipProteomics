#!/usr/bin/env Rscript

# title: supplment.R
# description: code associated with supplement.Rnw
# author: twab

## ---- fit0, results='hide'

## fit the protein-level model to WASHC4

# load dependencies
library(dplyr)
library(lmerTest)

# load SwipProteomics
data(swip)
data(msstats_prot)

# LMM formula
fx0 <- 'Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)'

# fit the model
fm0 <- lmer(fx0, msstats_prot %>% subset(Protein == swip))

# examine the model's summary
summary(fm0, ddf = "Satterthwaite")


## ---- contrast7 

## Compare 'Mutant:F7' and 'Control:F7' Conditions

# create a contrast
coeff <- lme4::fixef(fm0)
contrast7 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast7["GenotypeMutant:BioFractionF7"] <- +1 # positive coeff
contrast7["GenotypeControl:BioFractionF7"] <- -1 # negative coeff

# evaluate contrast
lmerTestContrast(fm0, contrast7)


## ---- contrast8 

# create a contrast to compare 'Mutant' versus 'Control'
contrast8 <- getContrast(fm0, "Mutant","Control")

# evaluate contrast
lmerTestContrast(fm0, contrast8)


## ---- fit1

# the module-level formula to be fit:
fx1 <- 'Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)'

# load WASH Complex proteins
data(washc_prots)

fm1 <- lmer(fx1, msstats_prot %>% subset(Protein %in% washc_prots))

# assess 'Mutant-Control' comparison
lmerTestContrast(fm1, contrast8)


## ---- nakagawa

# assess gof with Nakagawa coefficient of determination
r.squaredGLMM.merMod(fm0)

r.squaredGLMM.merMod(fm1)


## ---- variancePartition

# load variancePartition
library(variancePartition)

# calculate partitioned variance
form <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
fit <- lmer(form, data = msstats_prot %>% filter(Protein %in% washc_prots))

calcVarPart(fit)
