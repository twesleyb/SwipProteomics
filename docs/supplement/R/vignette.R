#!/usr/bin/env Rscript

# work though lmerTestContrast for protein- and module- level comparisions
# assess effect of setting intercept as +1 or +0
# assess Module-level fit with Protein as mixed or fixed effect

# * fit two types of models:
#     (1) protein-level model (fx1 -> fm1)
#     (2) module-level model (fx2 -> fm2)

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(swip) 
data(washc_prots)
data(msstats_prot)
data(msstats_contrasts)
data(mut_vs_control)

suppressPackageStartupMessages({
  library(lme4) # Bates et al., 2020
  library(lmerTest) # Kuznetsova et al. 2017
  library(dplyr)
  library(data.table)
  library(MSstatsTMT)
})

## 0. MSstatsTMT --------------------------------------------------------------

# intra-BioFraction comparisons
MSstatsTMT::groupComparisonTMT(msstats_prot %>% subset(Protein == swip),
			       contrast.matrix=msstats_contrasts) %>% knitr::kable()

##  lmerTestContrast = same result
fx <- Abundance ~ 1 + Condition + (1|Mixture)
fm <- lmer(fx, msstats_prot %>% subset(Protein == swip))

# create a contrast between coef(fm)
L8 <- getContrast(fm, "ConditionMutant.F8", "ConditionControl.F8")

# for a single intra-fraction comparison
lmerTestContrast(fm,L8) %>% knitr::kable()

# create a contrast for the Mutant-Control comparison
LT <- fixef(fm)
LT[] <- 0
LT[grep("Mutant",names(LT))] <- +1/7 # 7x
LT[grep("Control",names(LT))] <- -1/7 # 6x

# NOTE: the intercept is 0 and absorbs the first coeff==ConditionControl.F10
#length(names(which(LT < 0))) == 6 # Control.BioFraction
#length(names(which(LT > 0))) == 7 # Mutant.BioFraction
stopifnot(("ConditionControl.F10" %in% names(which(LT < 0)) ) == FALSE)
stopifnot(sum(LT)  == + 1/7)

lmerTestContrast(fm, LT) %>% mutate(Contrast="Mutant-Control") %>% 
	unique() %>% knitr::kable()

# overall Mut-Control comparison with MSstatsTMT
MSstatsTMT::groupComparisonTMT(msstats_prot %>% subset(Protein == swip),
			       contrast.matrix=mut_vs_control) %>% knitr::kable()

# NOTE: the DF differ subtly because of numerical precision term added to
# MSstatsTMT result.


## 1. fit protein-level model to WASHC4 ---------------------------------------

message("\nFitting protein-level mixed-models.")

fx <- "Abundance ~ Condition + (1|Mixture)"
fx0 <- "Abundance ~ 0 + Condition + (1|Mixture)"
fx1 <- "Abundance ~ 1 + Condition + (1|Mixture)"

# no intercept
fm <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein == swip))

summary(fm, ddf="Satterthwaite")

# intercept +0
message("\nfit: ",fx0)
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% subset(Protein == swip))

summary(fm0, ddf="Satterthwaite")

# intercept +1
message("\nfit: ",fx1)

fm1 <- lmerTest::lmer(fx1, msstats_prot %>% subset(Protein == swip))

summary(fm1, ddf="Satterthwaite")


## intra-biofraction contrasts

message("\nAssessing protein-level intra-BioFraction comparisons.")

# no intercept
message("\nfit: ",fx)

L7 <- getContrast(fm,"ConditionMutant.F7",
		  "ConditionControl.F7")

lmerTestContrast(fm1, L7) %>% mutate(Contrast="F7") %>% knitr::kable()

# +0
message("\nfit: ",fx0)

L7 <- getContrast(fm0,"ConditionMutant.F7",
		  "ConditionControl.F7")

lmerTestContrast(fm0,L7) %>% mutate(Contrast="F7") %>% knitr::kable()

# +1
message("\nfit: ",fx1)

L7 <- getContrast(fm1,"ConditionMutant.F7",
		  "ConditionControl.F7")

lmerTestContrast(fm1,L7) %>% mutate(Contrast="F7") %>% knitr::kable()


# illustrate the tricky comparison
# ConditionControl.F10 is embodied by the term (Intercept)
fm

# no intercept
L10 <- fixef(fm)
L10[] <- 0
L10["ConditionMutant.F10"] <- 1
#L10["ConditionControl.F10"] # doesnt exist
message("\nfit: ",fx)
lmerTestContrast(fm, L10) %>% mutate(Contrast="F10") %>% knitr::kable()

# easist with fm0
L10 <- getContrast(fm0, "ConditionMutant.F10","ConditionControl.F10")
message("\nfit: ", fx0)
lmerTestContrast(fm0, L10) %>% mutate(Contrast="F10") %>% knitr::kable()

# +1
L10 <- fixef(fm1)
L10[] <- 0
L10["ConditionMutant.F10"] <- 1
#L10["ConditionControl.F10"] # doesnt exist
message("\nfit:",fx1)
lmerTestContrast(fm1,L10) %>% mutate(Contrast="F10") %>% knitr::kable()


## 2. assess 'Mutant-Control' contrast -----------------------------------------

message("Assessing 'Mutant-Control' comparison.")

# no intercept
LT <- fixef(fm)
LT[] <- 0
LT[grep("Control",names(LT))] <- -1/7
LT[grep("Mutant",names(LT))] <- +1/7
message("\nfit: ",fx)
lmerTestContrast(fm, LT) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()

# easiest for fm0, create a contrast with getContrast
LT <- getContrast(fm0, "Mutant", "Control")
message("\nfit: ",fx0)
lmerTestContrast(fm0, LT) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()

# getContrast doesnt work for models with intercepts
LT <- fixef(fm1)
LT[] <- 0
LT[grep("Control",names(LT))] <- -1/7
LT[grep("Mutant",names(LT))] <- +1/7
message("\nfit: ",fx1)
lmerTestContrast(fm1, LT) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()


## - MSstatsTMT -------------------------------------------------------------------

# create mutant-Control Comparison
LT <- setNames(rep(0,ncol(msstats_contrasts)),colnames(msstats_contrasts))
LT[grepl("Mutant",names(LT))] <- +1/7
LT[grepl("Control",names(LT))] <- -1/7
msstats_contrasts <- rbind(msstats_contrasts,LT)
rownames(msstats_contrasts)[nrow(msstats_contrasts)] <- "Mutant-Control"

message("\nMSstatsTMT results:")
MSstatsTMT::groupComparisonTMT(msstats_prot %>% subset(Protein == swip), 
			       msstats_contrasts, moderated = FALSE) %>% 
  knitr::kable()


## 3. fit module-level model to the WASH complex -------------------------------

message("\nFitting module-level mixed-models.")

## the module-level model to be fit:
# three models with only subtly different results
fx <- "Abundance ~ Condition + Protein + (1|Mixture)" 
fx0 <- "Abundance ~ 0 + Genotype + BioFraction + Protein + (1|Mixture)"
fx1 <- "Abundance ~ 1 + Genotype + BioFraction + Protein + (1|Mixture)"

## Mutant-Control comparision
message("\nfit: ", fx)
fm <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein %in% washc_prots))
summary(fm,ddf="Satterthwaite")

message("\nfit: ", fx0)
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein %in% washc_prots))
summary(fm0,ddf="Satterthwaite")

message("\nfit: ", fx1)
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein %in% washc_prots))
summary(fm1,ddf="Satterthwaite")


## assess the overall comparison between Mutant and Control

# for fm, no intercept
LT <- fixef(fm)
LT[] <- 0
LT[grep("Control",names(LT))] <- -1/7
LT[grep("Mutant",names(LT))] <- +1/7
message("\nfit: ",fx)
lmerTestContrast(fm, LT) %>% mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# for +0
LT <- getContrast(fm0,"Mutant","Control") # stats differ
message("\nfit: ", fx0)
lmerTestContrast(fm0, LT) %>% mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

## getContrast doesnt work for models with +1 intercept
# for fm1: must specify contrast correctly!
LT <- fixef(fm1)
LT[] <- 0
LT[grep("Mutant",names(LT))] <- +1/7
LT[grep("Control",names(LT))] <- -1/7
message("\nfit: ", fx1)
lmerTestContrast(fm1, LT) %>% mutate(Contrast='Mutant-Control') %>% 
	unique() %>% knitr::kable()


## protein as a mixed-effect

fx_mixed <- "Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)"
fm_mixed <- lmerTest::lmer(fx_mixed, msstats_prot %>% filter(Protein %in% washc_prots))

fx_fixed <- "Abundance ~ 0 + Condition + Protein + (1|Mixture)"
fm_fixed <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein %in% washc_prots))

summary(fm_mixed)

summary(fm_fixed)

LT <- getContrast(fm_mixed,"Mutant","Control")
lmerTestContrast(fm_mixed, LT) %>% 
	mutate(Contrast='Mutant-Control') %>% unique() %>% 
	mutate(Protein = "mixed") %>% select(-isSingular) %>%
	knitr::kable()

LT <- getContrast(fm_fixed,"Mutant","Control")
lmerTestContrast(fm_fixed, LT) %>% 
	mutate(Contrast='Mutant-Control') %>% unique() %>% 
	mutate(Protein = "fixed") %>% select(-isSingular) %>%
	knitr::kable()

# numerically these results are the same


# try to compare the models
REMLcrit(fm_mixed) # protein as mixed

REMLcrit(fm_fixed) # protein as fixed

REMLcrit(fm_mixed) > REMLcrit(fm_fixed)

AIC(fm_mixed,fm_fixed) %>% knitr::kable()

AIC(fm_mixed) < AIC(fm_fixed)

r.squaredGLMM.merMod(fm_fixed)

# very similar, total R2c is greater for fm_mixed
r.squaredGLMM.merMod(fm_mixed)

# REML
# AIC
# nakagawa

# what about other modules?
data(partition)

module <- sample(split(names(partition),partition),1)
M <- paste0("M",names(module))
prots <- unlist(module,use.names=FALSE)

cbind(M,nProts=length(prots))

fx_mixed <- "Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)"
fm_mixed <- lmerTest::lmer(fx_mixed, msstats_prot %>% filter(Protein %in% prots))
fx_fixed <- "Abundance ~ 0 + Condition + Protein + (1|Mixture)"
fm_fixed <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein %in% prots))

# M61

# fixed
r.squaredGLMM.merMod(fm_fixed) %>% knitr::kable() # nearly all variance attributable to fixef

REMLcrit(fm_fixed)

AIC(fm_fixed)

# mixed
r.squaredGLMM.merMod(fm_mixed) %>% knitr::kable() 
# fixef now only account for 48% of total -- some of this variance is now
# attributable to mixed effects. modeling protein as a mixed effect we can
# seperate variance attributable to Condition and Protein

REMLcrit(fm_mixed)

AIC(fm_mixed)
