#!/usr/bin/env Rscript

# 1. formula("Abundance ~ 0 + Condition + (1 | Mixture)")
# 2. formula("Abundance ~ 0 + Genotype + BioFraction + (1 | Subject)")
# 3. formula("Abundance ~ 0 + Genotype + ... + (1|Protein)")

## prepare the environment
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)


## load project data
data(swip)
data(fx0); data(cm0)
data(fx1); data(cm1)
data(msstats_prot)
data(msstats_contrasts)


## imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MSstatsTMT)
})


## illustrate analysis with lmerTest

# fit a model
fm <- lmerTest::lmer(fx1,data=msstats_prot %>% filter(Protein == swip))

# given the model, test a single contrast:
lmerTestContrast(fm, cm1) %>% knitr::kable()

## or assess multiple contrasts with lmerTestProtein
# FIXME: need to add capability for multiple proteins 
# FIXME: need to add moderated capability
lmerTestProtein(swip, fx0, msstats_prot, cm0) %>% knitr::kable() 

# NOTE: cm0 is a list of numeric vectors indicating positive and negative coeff
# defining a comparison to be tested

# we can also pass msstats_contrasts (a matrix) to lmerTestProtein:
lmerTestProtein(swip, fx0, msstats_prot, msstats_contrasts) %>% knitr::kable() 


## [1] fx0 / fitBioFraction / intra-fraction comparisons -----------------------

# store results in a list
results <- list()

## analysis with lmerTest
input = list(swip, fx0, msstats_prot, cm0) 
results[["lmerTest::IntraFraction"]] <- do.call(lmerTestProtein, input)

## analysis with MSstats 
input = list("data" = msstats_prot %>% filter(Protein==swip),
	     "contrast.matrix" = msstats_contrasts,
	     "moderated"=FALSE)
results[["MSstatsTMT::IntraFraction"]] <- suppressMessages({
	do.call(groupComparisonTMT, input) })

# FIXME: Remove message printing formula


# [2] fx1 -- fitGenotype - Mutant vs Control ----------------------------------

## analysis with lmerTest
input = list(swip, fx1, msstats_prot, cm1) 
results[["lmerTest::Mutant-Control"]] <- do.call(lmerTestProtein, input)

## analysis with MSstats 

# create a contrast for assessing difference between Control and Mutant
alt_contrast <- matrix(c(-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,
		     1/7,1/7,1/7,1/7,1/7,1/7,1/7), nrow=1)
row.names(alt_contrast) <- "Mutant-Control"
colnames(alt_contrast)<- levels(msstats_prot$Condition)


# do MSstatsTMT groupComparisons
input = list("data" = msstats_prot %>% filter(Protein==swip),
	     "contrast.matrix" = alt_contrast,
	     "moderated"=FALSE)
results[["MSstatsTMT::Mutant-Control"]] <- suppressMessages({
	do.call(groupComparisonTMT, input) })


###############################################################################

# NOTE: the above approaches differ in the way the model and contrasts are
# specified. Will lmerTest approach return the same results?

###############################################################################


## Control-Mutant comparisons with alternative contrast and lmerTest ----------
# the result is not the same. Something is different here.
# somehow the result is different, even though we expected the same result

## analysis with lmerTest
input = list(swip, fx0, msstats_prot, alt_contrast) 
results[["lmerTest::alt-Mutant-Control"]] <- do.call(lmerTestProtein, input)

# ^this result matches MSstats

## Examine results:
lapply(results,knitr::kable)


###############################################################################

# which is correct way to specify the model?

# as ther are no sig prots with lmerTestContrast(fm1), it seems like something
# may be wrong

###############################################################################



## timed comparison -----------------------------------------------------------

# do MSstatsTMT groupComparisons
MSstatsTMT <- function() { 
	input = list("data" = msstats_prot %>% filter(Protein==swip),
		     "contrast.matrix" = msstats_contrasts,"moderated"=FALSE)
        results <- suppressMessages({ do.call(groupComparisonTMT, input) })
}

# do lmerTestProtein
lmerTest <- function() {
	input = list(swip, fx0, msstats_prot, msstats_contrasts) 
	results <- do.call(lmerTestProtein, input)
}

timed_res <- microbenchmark(MSstatsTMT(), lmerTest(), times=100L)   

print(timed_res)

# it appears that lmerTest approach is significantly faster ~1.5x
df <- as.data.table(do.call(cbind,timed_res)) 

# time to evaluate 10,000 proteins
fold_diff <- df %>% group_by(expr) %>% 
	summarize(mean=mean(time),.groups="drop") %>%
	mutate("Duration" = 10000*(mean*10^-9)/60)

# savings of ~7 min
