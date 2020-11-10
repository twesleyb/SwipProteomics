#!/usr/bin/env Rscript

#library(merDeriv)
library(SwipProteomics)

devtools::load_all()

data(swip)
data(gene_map)
data(msstats_prot)

## ----------------------------------------------------------------------------

## the protein-level model to be fit:
fx0 <- formula(Abundance ~ 0 + Genotype:BioFraction + (1|Mixture))
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

l1 <- getContrast(fm0, "Control", "Mutant")

results <- lmerTestContrast(fm0, l1, variance=TRUE) 
results %>% knitr::kable()

attr(results,"variance") %>% knitr::kable()


## ----------------------------------------------------------------------------
## the protein-level model to be fit:
fx0 <- formula(Abundance ~ 0 + Genotype:BioFraction + (1|Mixture))
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

# calculate variance partition
#FIXME: need to explain why modeled all as mixef
form <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture)"
fit <- lme4::lmer(form, msstats_prot %>% filter(Protein == swip))
mixef_var <- as.data.frame(lme4::VarCorr(fit,comp="Variance"))
x <- setNames(mixef_var$vcov,nm=mixef_var$grp) # variance of each component
PVE_swip = x/sum(x)

#sum(PVE_swip)==1

## ----------------------------------------------------------------------------

## the module-level model to be fit:
washc_prots = mapID("Washc*")
form <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"

# warnings because our module is too perfect
fit <- suppressWarnings({ lmerTest::lmer(form, msstats_prot %>% filter(Protein %in% washc_prots)) })

# calculate partitioned variance
mixef_var <- as.data.frame(lme4::VarCorr(fit,comp="Variance"))
x <- setNames(mixef_var$vcov,nm=mixef_var$grp) # variance of each component
PVE_washc = x/sum(x)
PVE_washc %>% knitr::kable()
