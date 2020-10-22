#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(partition)
data(msstats_prot)
data(fx0)
data(fx1)

library(dplyr)
library(data.table)

msstats_prot$Module <- paste0("M",partition[msstats_prot$Protein])

msstats_filt <- msstats_prot %>% filter(!is.na(Module))


fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Protein)")

module = "M1"

subdat <- msstats_filt %>% filter(Module==module)

fm <- lmerTest::lmer(fx,subdat)

summary(fm,ddf="Satterthwaite")

r.squaredGLMM.merMod(fm) %>% knitr::kable()
