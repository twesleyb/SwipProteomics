#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)
devtools::load_all(root,quiet=TRUE)

data(gene_map)
data(ne_surprise_partition)
data(msstats_prot)

library(dplyr)
library(data.table)

library(geneLists)

data(corum)

corum["Fibrinogen complex"]

entrez <- corum[["Fibrinogen complex"]]
fibrinogen <- mapID(entrez,"entrez","uniprot")

save(fibrinogen,file=file.path(root,"data","fibrinogen.rda"),version=2)


fx = Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein)

fm = lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% fibrinogen))
summary(fm)

LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT) %>% mutate(Contrast="Mutant-Control") %>% unique() %>% knitr::kable()
