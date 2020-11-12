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
lmerTestContrast(fit2, L2) %>% mutate(Contrast = 'Mutant-Control') %>% unique() %>% knitr::kable()

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

###################################################################################

# try glmer

library(lattice)

xyplot(incidence/size ~ period|herd, cbpp, type=c('g','p','l'),layout=c(3,5), index.cond = function(x,y)max(y))

library(lme4)

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial())


data(msstats_prot)
data(washc_prots)
fx <- "scale_Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)"
subdat <- msstats_prot %>% group_by(Protein) %>% mutate(scale_Abundance = Abundance/max(Abundance)) %>% ungroup() %>% subset(Protein %in% washc_prots)
gm1 <- glmer(fx, data = subdat, family = poisson())

fx <- "Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)"
subdat <- msstats_prot %>% group_by(Protein) %>% mutate(scale_Abundance = Abundance/max(Abundance)) %>% ungroup() %>% subset(Protein %in% washc_prots)
gm1 <- glmer(fx, data = subdat, family = poisson())


