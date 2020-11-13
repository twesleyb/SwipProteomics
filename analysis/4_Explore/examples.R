#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)

#library(merDeriv)
library(dplyr)
#library(SwipProteomics)

devtools::load_all()

data(swip)
data(gene_map)
data(msstats_prot)

## ----------------------------------------------------------------------------

## the module-level model to be fit:
washc_prots = mapID("Washc*")

# variance partitioned
form0 <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
fit0 <- lmerTest::lmer(form0, msstats_prot %>% filter(Protein %in% washc_prots))
vp <- variancePartition::calcVarPart(fit0)
sum(vp) == 1

## these results are the same!

# calculate partitioned variance
var_df <- as.data.frame(lme4::VarCorr(fit0,comp="Variance"))
mixef_var <- setNames(var_df$vcov,nm=var_df$grp) # variance of each component
round(mixef_var/sum(mixef_var),4)


## Mutant-Control comparision
# two models with the same result
form2 <- "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)"
fit2 <- lmerTest::lmer(form2, msstats_prot %>% filter(Protein %in% washc_prots))
L2 <- getContrast(fit2,"Mutant","Control")
lmerTestContrast(fit2, L2) %>% mutate(Contrast = 'Mutant-Control') %>% unique() %>% knitr::kable()

form1 <- "Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)"
fit1 <- lmerTest::lmer(form1, msstats_prot %>% filter(Protein %in% washc_prots))
L1 <- getContrast(fit1,"Mutant","Control")
lmerTestContrast(fit1,L1)  %>% knitr::kable()

getVariance(fit1) # sigma^2


qqnorm(resid(fit0))
qqline(resid(fit0))

qqnorm(resid(fit1))
qqline(resid(fit1))

qqnorm(resid(fit2))
qqline(resid(fit2))


## real world module

data(partition)

modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

wash_module <- paste0("M",partition[swip])
prots <- modules[[wash_module]]

# variance partitioned
fit0 <- lmerTest::lmer(form0, msstats_prot %>% filter(Protein %in% prots))
vp <- variancePartition::calcVarPart(fit0)
vp # 81 % of variance is attributable to protein
# 5% of variance is attributable to Genotype

## these results are the same!

# calculate partitioned variance
var_df <- as.data.frame(lme4::VarCorr(fit0,comp="Variance"))
mixef_var <- setNames(var_df$vcov,nm=var_df$grp) # variance of each component
round(mixef_var/sum(mixef_var),4)


## Mutant-Control comparision

# two models with the same result
form2 <- "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)"
fit2 <- lmerTest::lmer(form2, msstats_prot %>% filter(Protein %in% prots))
L2 <- getContrast(fit2,"Mutant","Control")
lmerTestContrast(fit2, L2) %>% mutate(Contrast = 'Mutant-Control') %>% 
	unique() %>% knitr::kable()


form1 <- "Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)"
fit1 <- lmerTest::lmer(form1, msstats_prot %>% filter(Protein %in% prots))
L1 <- getContrast(fit1,"Mutant","Control")
lmerTestContrast(fit1,L1)  %>% knitr::kable()

## fat tails

qqnorm(resid(fit0))
qqline(resid(fit0))

qqnorm(resid(fit1))
qqline(resid(fit1))

qqnorm(resid(fit2))
qqline(resid(fit2))


## fit to proteasome
library(geneLists)

data(corum)
data(gene_map)
data(partition)
data(msstats_prot)

idx <- which(names(corum) == "26S proteasome")
proteasome <- mapID(corum[[idx]],"entrez","uniprot")

sum(proteasome %in% names(partition))

fx1 <- "Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)"
fm <- lmerTest::lmer(fx1, msstats_prot %>% subset(Protein %in% proteasome)) 

var1 <- getVariance(fm)
R2 <- var1["fixed"]/sum(var1) # nakagawa!

LT <- getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>% knitr::kable()

r.squaredGLMM.merMod(fm) %>% knitr::kable()

module <- names(which(partition==83))
fm83 <- lmerTest::lmer(fx1, data = msstats_prot %>% subset(Protein %in% module)) 

lmerTestContrast(fm83,LT) %>% knitr::kable()

# We need to know the PVE for each term!

library(variancePartition)

fx2 <- "Abundance ~  (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"

fm <- lmerTest::lmer(fx2, data = msstats_prot %>% subset(Protein %in% proteasome)) 
vp1 <- calcVarPart(fm)

# This module is highly cohesive. Protein variation only explains ~7% of
# variance.

qqnorm(resid(fm))
qqline(resid(fm))

# this is common... fat tails

module <- names(which(partition==83))
fm <- lmerTest::lmer(fx2, data = msstats_prot %>% subset(Protein %in% module)) 
vp2 <- calcVarPart(fm)

getVariance(fm)


# no variance explained by genotype
x <- as.numeric(vp1["BioFraction"] + vp1["Genotype"])
y <- as.numeric(vp1["Protein"])

x # maximize
y # minimize

q1 <- x/y

x <- as.numeric(vp2["BioFraction"] + vp2["Genotype"])
y <- as.numeric(vp2["Protein"])
q2 <- x/y



data(washc_prots)

fx <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
fm <- lmerTest::lmer(fx, data = msstats_prot %>% subset(Protein %in% washc_prots))
vp0 <- calcVarPart(fm)
	
partition <- setNames(c(1,1,1,2,2),nm=washc_prots)
modules <- split(partition,partition)

fit_list <- list()
for (i in seq(modules)){
	fit_list[[i]] <- lmerTest::lmer(fx, data = msstats_prot %>% subset(Protein %in% names(modules[[i]])))
}

q0 <- as.numeric((vp0["Genotype"] + vp0["BioFraction"]) / vp0["Protein"])

vp1 <- calcVarPart(fit_list[[1]])
q1 <- as.numeric((vp1["Genotype"] + vp1["BioFraction"]) / vp1["Protein"])

vp2 <- calcVarPart(fit_list[[2]])
q2 <- as.numeric((vp2["Genotype"] + vp2["BioFraction"]) / vp2["Protein"])


