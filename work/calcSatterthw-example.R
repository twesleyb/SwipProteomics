#!/usr/bin/env Rscript

# title: calcSatterth-example.R
# description: usage of calcSattert

## ---- prepare the renv

root = "~/projects/SwipProteomics"
renv::load(root)

devtools::load_all(root)


## ---- functions

cleanSatterthwaiteTable <- function(data, ...) {
	data %>% unlist() %>% 
		t() %>% 
		as.data.table() %>% 
		mutate(pvalue = formatC(pvalue)) %>%
		knitr::kable()
}

cleanContrastTable <- function(data, ...) {
	data %>% mutate(Contrast = "Mutant-Control") %>%
		mutate(Pvalue = formatC(Pvalue)) %>%
		select(-isSingular) %>%
		unique() %>% knitr::kable()
}


## ---- imports

suppressPackageStartupMessages({
  library(dplyr)
  library(lmerTest)
  library(data.table)
})


## ---- load the data

# library(SwipProteomics)
data(washc_prots)
data(msstats_prot)


## ---- fit model to WASH Complex

# the lmer function with mixef of Mixture and Protein
fx <- Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein)

# fit the model
fm <- lmer(fx, msstats_prot %>% filter(Protein %in% washc_prots))

summary(fm)


## ---- create a contrast

# a comparison between Mutant and Control coefficients
LT <- getContrast(fm, "Mutant","Control")


## ---- lmerTestContrast

lmerTestContrast(fm, LT) %>% cleanContrastTable()


## ---- lmerTest::calcSatterth

# the degrees of freedome and Pvalue are the same
# the F-statistic is reported instead of the t-statistic
lmerTest::calcSatterth(fm,LT) %>% cleanSatterthwaiteTable()
