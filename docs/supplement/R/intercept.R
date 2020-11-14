#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: set the intercept to +1

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

fx0 <- "Abundance ~ 0 + Condition + (1|Mixture)"
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% subset(Protein == swip))

df <- summary(fm0, ddf="Satterthwaite")[["coefficients"]]
df <- as.data.table(df,keep.rownames="Coefficient")
colnames(df)[colnames(df) == "Pr(>|t|)"] <- "p value"
df$"p value" <- formatC(df$"p value")
df %>% knitr::kable()

L8 <- getContrast(fm0,"Mutant","Control")
lmerTestContrast(fm0,L8) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()

# fit protein-level model
# NOTE: we specify the intercept = 1
fx <- "Abundance ~ 1 + Condition + (1|Mixture)"
fm <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein == swip))

# examine coefficients
df <- summary(fm, ddf="Satterthwaite")[["coefficients"]]
df <- as.data.table(df,keep.rownames="Coefficient")
colnames(df)[colnames(df) == "Pr(>|t|)"] <- "p value"
df$"p value" <- formatC(df$"p value")
df %>% knitr::kable()

# I find it much harder to interpret this table...
L8 <- fixef(fm)
L8[] <- 0
L8[which(grepl("Mutant",names(L8)))] <- +1/7
L8[which(grepl("Control",names(L8)))] <- -1/7

lmerTestContrast(fm,L8) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()

quit()

# we must specify the contrast correctly.
L10 <- fixef(fm)
L10[] <- 0

# The tricky one is the contrast for ConditionMutant.F10 versus
# ConditionControl.F10. --> the intercept represents ConditionControl.F10?

L10["ConditionMutant.F10"] <- 1
L10["ConditionControl.F10"] # doesnt exist
L10

lmerTestContrast(fm,lT) %>% knitr::kable() # this is correct

# try one more
L4 <- L10
L4[] <- 0
L4["ConditionMutant.F4"] <- 1
L4["ConditionControl.F4"] <- -1 

lmerTestContrast(fm,L4) %>% knitr::kable() # this is correct


## fit to wash module ---------------------------------------------------------

data(partition)
data(washc_prots)


# NOTE: we specify the intercept = 1
fx1 <- "Abundance ~ 1 + Condition + Protein + (1|Mixture)"
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% subset(Protein %in% washc_prots))

qqnorm(residuals(fm1))
qqline(residuals(fm1))


# specify the overall contrast.
LT <- fixef(fm1)
LT[] <- 0
LT[grepl("Control",names(LT))] <- -1/7
LT[grepl("Mutant",names(LT))] <- +1/7

lmerTestContrast(fm1,LT) %>% knitr::kable()


## expand scope

modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

#partition[swip] # "M17"
m = "M17"
prots <- modules[[m]]
fit <- lmerTest::lmer(fx1, msstats_prot %>% subset(Protein %in% prots))

# specify the overall contrast.
LT <- fixef(fit)
LT[] <- 0
LT[grepl("Control",names(LT))] <- -1/7
LT[grepl("Mutant",names(LT))] <- +1/7

lmerTestContrast(fit,LT) %>% knitr::kable()

qqnorm(residuals(fit))
qqline(residuals(fit))
