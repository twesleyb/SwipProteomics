#!/usr/bin/env Rscript

# title: supplment.R
# description: code associated with supplement.Rnw
# author: twab

## ---- fit0, results='hide'

# load dependencies
library(dplyr)
library(lmerTest)

# load SwipProteomics
data(swip)
data(msstats_prot)

# formula to be fit to WASHC4, aka SWIP:
fx0 <- 'Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)'

# fit the LMM
fm0 <- lmer(fx0, msstats_prot %>% subset(Protein == swip))

# examine the model's summary
summary(fm0, ddf = "Satterthwaite")


## ---- contrast0 
#<<contrast0, eval=T, echo=T, results='hide'>>=

## Compare 'Mutant:F7' and 'Control:F7'

# create a contrast
coeff <- lme4::fixef(fm0)
contrast7 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast7["GenotypeMutant:BioFractionF7"] <- +1 # positive coeff
contrast7["GenotypeControl:BioFractionF7"] <- -1 # negative coeff

# evaluate contrast
lmerTestContrast(fm0, contrast7)


## Compare 'Mutant' versus 'Control'

# create a contrast
contrast8 <- getContrast(fm0, "Mutant","Control")

# evaluate contrast
lmerTestContrast(fm0, contrast8)


## ---- module-model
# <<module-model, eval=TRUE, echo=TRUE, results='hide'>>=

# the module-level formula to be fit:
fx1 <- 'Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)'

# load WASH Complex proteins
data(washc_prots)

fm1 <- lmer(fx1, msstats_prot %>% subset(Protein %in% washc_prots))


## ---- rsquared
#<<rsquared, eval=TRUE, echo=TRUE, results='hide'>>=

# assess gof with Nakagawa coefficient of determination
r.squaredGLMM.merMod(fm1)

r.squaredGLMM.merMod(fm0)


## ---- variancePartition
#<<variancePartition, eval=TRUE, echo=FALSE, results='hide'>>=

# load variancePartition
library(variancePartition)

# calculate partitioned variance
form <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
fit <- lmer(form, data = msstats_prot %>% filter(Protein %in% washc_prots))

calcVarPart(fit)
