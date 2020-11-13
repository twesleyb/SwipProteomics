#!/usr/bin/env Rscript

# work though lmerTestContrast for protein- and module- level comparisions

# * fit two types of models:
#     (1) protein-level model (fx1 -> fm1)
#     (2) module-level model (fx2 -> fm2)

# * analyze the variance partition of each model 

root <- "~/projects/SwipProteomics"
renv::load(root)


library(lme4) # Bates et al., 2020


library(lmerTest) # Kuznetsova et al. 2017


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

#renv::install("twesleyb/SwipProteomics") # github
#library(SwipProteomics)
devtools::load_all(root)

data(swip) 
data(msstats_prot)

# examine key elements of msstats_prot
msstats_prot %>% select(Abundance, Protein, Condition) %>% head() %>% 
	knitr::kable()


## 1. fit protein-level model to WASHC4 ---------------------------------------

# fit protein-level model
fx1 <- "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)"
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% subset(Protein == swip))

(fm1)

# examine coefficients
df <- summary(fm1,ddf="Satterthwaite")[["coefficients"]]
df <- as.data.table(df,keep.rownames="Coefficient")
colnames(df)[colnames(df) == "Pr(>|t|)"] <- "p value"
df$"p value" <- formatC(df$"p value")
df %>% knitr::kable()

L7 <- getContrast(fm1,"GenotypeMutant:BioFractionF7",
		  "GenotypeControl:BioFractionF7")

lmerTestContrast(fm1,L7) %>% mutate(Contrast="F7") %>% knitr::kable()


## 2. assess 'Mutant-Control' contrast -----------------------------------------


L8 <- getContrast(fm1,"Mutant","Control")


results <- lmerTestContrast(fm1, L8) 


# examine results
results %>% select(-isSingular) %>% 
	mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()


## 3. fit module-level model to the WASH complex -------------------------------

data(gene_map)

## the module-level model to be fit:
washc_prots <- mapID("Washc*")

## Mutant-Control comparision

# two models with the same result
fx2 <- "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)"

fm2 <- lmerTest::lmer(fx2, msstats_prot %>% filter(Protein %in% washc_prots))

# examine the model
summary(fm2, ddf="Satterthwaite")

## assess the overall comparison between Mutant and Control
L8 <- getContrast(fm2,"Mutant","Control")

results <- lmerTestContrast(fm2, L8)

# examine the results
results %>% select(-isSingular) %>% 
	mutate(Contrast = 'Mutant-Control') %>% unique() %>% knitr::kable()


## 4. alternative model -------------------------------------------------------
# alternative model with the same result

alt_fx2 <- "Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)"
alt_fm2 <- lmer(alt_fx2, msstats_prot %>% filter(Protein %in% washc_prots))

L8 <- getContrast(alt_fm2, "Mutant","Control")

results <- lmerTestContrast(alt_fm2, L8)  

results %>% select(-isSingular) %>% mutate(Contrast="Mutant-Control") %>% 
	unique() %>% knitr::kable()



## 5. variancePartition -------------------------------------------------------

# NOTE: fixed effects account for most of the variation
getVariance(alt_fm2) # sigma^2

# library(variancePartition)

# variance partitioned
vp_fx <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
vp_fm <- lmerTest::lmer(vp_fx, msstats_prot %>% filter(Protein %in% washc_prots))

vp <- variancePartition::calcVarPart(vp_fm)

knitr::kable(t(vp))
sum(vp) == 1

## these results are the same!

# calculate partitioned variance
var_df <- VarCorr(vp_fm, comp="Variance") %>% as.data.table()
mixef_var <- setNames(var_df$vcov, nm=var_df$grp) # variance of each component

round(mixef_var/sum(mixef_var),4) %>% t() %>% knitr::kable()


## qqplot protein-model fm1 ---------------------------------------------------

# fit to WASHC4
qqnorm(resid(fm1))
qqline(resid(fm1))


## qqplot module-model fm2 ---------------------------------------------------

# fit to WASH Complex
qqnorm(resid(fm2))
qqline(resid(fm2))

# fit to WASH Complex
qqnorm(resid(alt_fm2))
qqline(resid(alt_fm2))


## 6. real world module -------------------------------------------------------

# module  that contains the wash complex protein, SWIP

data(partition) 

# vector of module membership for all proteins
head(partition)

# all modules
modules <- split(partition,partition)

# proteins in same module as WASHC4 (swip)
prots <- names(which(partition==partition[swip]))

length(prots)


## Mutant-Control comparision

# alternative model, same result
fx2 <- "Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)"

fit <- lmerTest::lmer(fx2, msstats_prot %>% filter(Protein %in% prots))

L8 <- getContrast(fit,"Mutant","Control")

lmerTestContrast(fit,L8) %>% select(-isSingular) %>% 
	mutate(Contrast='Mutant-Control') %>% knitr::kable()


## proteasome example module ---------------------------------------------------

# load SwipProteomics data
data(gene_map)
data(partition)
data(msstats_prot)

# get proteomse gene list
library(geneLists) # twesleyb/geneLists

# load corum gene lists (mapped to mouse entrez)
data(corum)


## fit to proteasome

idx <- which(names(corum) == "26S proteasome")
proteasome <- mapID(corum[[idx]],"entrez","uniprot")

head(proteasome)

knitr::kable(cbind(module="proteasome",
		   nProts=sum(proteasome %in% names(partition))))

# fit the alternative module-level model
fit <- lmerTest::lmer(alt_fx2, msstats_prot %>% subset(Protein %in% proteasome)) 

(fit)

# compute sigma
var1 <- getVariance(fit)

var1

R2 <- var1["Fixed"]/sum(var1) # nakagawa!

R2 # the total variance explained by fixed effects


## assess Mutant-Control contrast

L8 <- getContrast(fit,"Mutant","Control")

lmerTestContrast(fit, L8) %>% select(-isSingular) %>% 
	mutate(Contrast='Mutant-Control') %>% knitr::kable()

# no change

knitr::kable(t(var1)) # variance of mixed and fixed effects

# highly cohesive, variance explained by Protein is small.
# variance of Fixed effects, Condition (Genotype:BioFraction)
# accounts for a large proportion of variation within the module.


## variancePartition ----------------------------------------------------------

# NOTE: We need to know the PVE for each term!

#library(variancePartition)
vp_fx <- "Abundance ~  (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
vp_fm <- lmerTest::lmer(vp_fx, msstats_prot %>% subset(Protein %in% proteasome))

x <- getVariance(vp_fm) # The variance attributed to fixed effect terms is 0!
knitr::kable(t(x/sum(x))) # as a percentage of the total

# Genotype accounts for very little variation. BioFraction is the major source
# of variation for the fixed effects.  

# Interpretation:
# This module is highly cohesive. Protein variation only explains ~7% of
# variance. Variation from BioFraction dominates.

## goodness of fit to proteasome ----------------------------------------------

qqnorm(resid(fit))
qqline(resid(fit))
