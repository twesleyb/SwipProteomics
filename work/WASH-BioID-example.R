#!/usr/bin/env Rscript

## ---- renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)
devtools::load_all(root,quiet=TRUE)

data(swip)
data(bioid_prot)


## ---- partition variance

# check variance of covariates
fx <- Abundance ~ (1|Run) + (1|Condition) + (1|Subject)
fm <- lme4::lmer(fx,bioid_prot %>% subset(Protein == swip))
vp <- getVariance(fm)
pve <- vp/sum(vp)
pve

## ---- lmTestContrast

# simple linear model, int=0 (explicitly estimate all coeff)
fx <- Abundance ~ 0 + Condition
fm <- lm(fx,bioid_prot %>% subset(Protein == swip))

LT = getContrast(fm,"WASH","Control")
lmTestContrast(fm,LT)
