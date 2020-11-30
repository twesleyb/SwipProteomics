#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load_all(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)

data(gene_map)
data(partition)
data(msstats_prot)

glud1 <- mapID("Glud1","symbol","uniprot")
glud1

fx = Abundance ~ 1 + Condition + (1|Mixture)

fm = lmerTest::lmer(fx, msstats_prot %>% filter(Protein == glud1))
fm

LT = getContrast(fm, "Mutant","Control")
lmerTestContrast(fm,LT)
