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


## 1. fit protein-level model to WASHC4 ---------------------------------------

# fit protein-level model
fx1 <- "Abundance ~ Condition + (1|Mixture)"
fx2 <- "Abundance ~ 1 + Condition + (1|Mixture)"
fx3 <- "Abundance ~ 1 + Genotype + BioFraction + (1|Mixture)"
fx4 <- "Abundance ~ 1 + Genotype:BioFraction + (1|Mixture)" # rank deficient

args <- list(fx2, msstats_prot %>% subset(Protein == swip))
fm <- do.call(lmerTest::lmer,args)

fm1 <- lmerTest::lmer(fx2, msstats_prot %>% subset(Protein == swip))


# examine coefficients
df <- summary(fm1,ddf="Satterthwaite")[["coefficients"]]
df <- as.data.table(df,keep.rownames="Coefficient")
colnames(df)[colnames(df) == "Pr(>|t|)"] <- "p value"
df$"p value" <- formatC(df$"p value")
df %>% knitr::kable()

# hard to interpret...
