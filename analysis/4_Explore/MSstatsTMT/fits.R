#!/usr/bin/env Rscript

# 1. formula("Abundance ~ 0 + Condition + (1 | Mixture)")
# 2. formula("Abundance ~ 0 + Genotype + BioFraction + (1 | Subject)")
# 3. formula("Abundance ~ 0 + Genotype + ... + (1|Protein)")

# prepare the environment
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load project data
data(swip)
data(fx0); data(cm0)
data(fx1); data(cm1)
data(msstats_prot)
data(msstats_contrasts)

# imports
library(dplyr)


# [1] fx0 / fitBioFraction / intra-fraction comparisons -----------------------

## analysis with lmerTest

# fit a model
fm0 <- lmerTest::lmer(fx0,data=msstats_prot %>% filter(Protein == swip))

# you can assess a single contrast with lmerTestContrast
#lmerTestContrast(fm0, cm0[[1]]) %>% knitr::kable()

# use lmerTestProtein to assess all contrasts in a list or matrix.
# FIXME: need to add capability for multiple proteins 
# FIXME: need to add moderated capability
res0 <- lmerTestProtein(swip, fx0, msstats_prot, cm0) 

## lmerTest results
res0 %>% knitr::kable()

## analysis with MSstats 
msstats0 <- MSstatsTMT::groupComparisonTMT(msstats_prot %>% filter(Protein==swip),
			       contrast.matrix = msstats_contrasts,
			       moderated=FALSE)

msstats0 %>% knitr::kable()


# the results appear to be the same.

# we can also pass msstats_contrasts to lmerTestProtein:
foobar <- lmerTestProtein(swip, fx0, msstats_prot, msstats_contrasts) 

foobar %>% knitr::kable()


# [2] fx1 -- fitGenotype - Mutant vs Control ----------------------------------

## analysis with lmerTest
res1 <- lmerTestProtein(swip, fx1, msstats_prot, cm1) 
res1 %>% knitr::kable()

## analysis with MSstats 

# create a contrast
contrast <- matrix(c(-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,
		     1/7,1/7,1/7,1/7,1/7,1/7,1/7), nrow=1)
row.names(contrast) <- "Mutant-Control"
colnames(contrast)<- levels(msstats_prot$Condition)

msstats1 <- MSstatsTMT::groupComparisonTMT(msstats_prot %>% filter(Protein==swip),
			       contrast.matrix = contrast,
			       moderated=FALSE)

msstats1 %>% knitr::kable()

# NOTE: the above approaches differ in the way the model and contrasts are
# specified. Will lmerTest approach return the same results?

# first, create a contrast 
contrast <- lme4::fixef(fm0)
contrast[grep("Control",names(contrast))] <- -1/7
contrast[grep("Mutant",names(contrast))] <- +1/7

# Use fm0 and new contrast
res1_alt = lmerTestContrast(fm0,contrast) 

res1_alt %>% knitr::kable()

# the result is not the same. Something is different here.
# somehow the result is different, even though we expected the same result


#------------------------------------------------------------------------------


