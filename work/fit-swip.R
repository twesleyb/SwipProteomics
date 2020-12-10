#!/usr/bin/env Rscript

## ---- prepare the renv

root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)

library(dplyr)
library(data.table)


## ---- load the data

data(swip)
data(swip_tmt)
data(washc_prots)

# rel_Intensity data
swip_tmt <-  swip_tmt %>% 
	group_by(Protein) %>% 
	dplyr::mutate(rel_Intensity = Intensity/max(Intensity))


## ---- protein-level comparison 

fx0 <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)

# fit the model
message(fx0)
fm0 <- lmerTest::lmer(fx0, swip_tmt %>% subset(Protein == swip))
summary(fm0)

# get a contrast:
LT <- getContrast(fm0,"Mutant","Control")

# assess contrast
lmerTestContrast(fm0, LT) %>% 
	dplyr::mutate(Contrast='Mutant-Control') %>% 
	unique() %>% 
	knitr::kable()


## ---- Module-level fit and contrast

fx1 <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

# fit the model
message(fx1)
fm1 <- lmerTest::lmer(fx1, swip_tmt %>% subset(Protein %in% washc_prots))
summary(fm1)

# get a contrast:
LT <- getContrast(fm1,"Mutant","Control")

# assess contrast
lmerTestContrast(fm1, LT) %>% 
	dplyr::mutate(Contrast='Mutant-Control') %>% 
	unique() %>% 
	knitr::kable()
