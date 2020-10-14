#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root)

devtools::load_all()

data(swip)
data(msstats_prot)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(variancePartition)
})

data(samples)
data(msstats_prot)

# munge to split Condition (Geno.BioFrac) into geno annotation
samples$Genotype <- sapply(strsplit(samples$Condition,"\\."),"[",1)

# munge to create mouse Subject identifier
samples$Subject <- as.numeric(interaction(samples$Genotype,samples$Experiment))

# munge to annotate msstats_prot with Subject (1-6)
msstats_prot$Subject <- as.numeric(interaction(msstats_prot$Genotype,msstats_prot$Mixture))

# cast the data into a matrix.
fx <- formula(Protein ~ Mixture + Channel + BioFraction + Genotype + Subject)
dm <- msstats_prot %>%
	reshape2::dcast(fx, value.var= "Abundance") %>% 
	as.data.table() %>% na.omit() %>% as.data.table() %>% as.matrix(rownames="Protein")

# munge to create sample info from the dcast fx
info <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(info) <- strsplit(as.character(fx)[3]," \\+ ")[[1]]

# Create a formula specifying all covariates as mixed effects and do
# variancePartition
form = formula(~ (1|Mixture) + (1|Channel) + (1|BioFraction) + (1|Genotype) + (1|Subject))
prot_varpart = fitExtractVarPartModel(dm, form, info)

# save
myfile <- file.path(root,"rdata","prot_varpart.rda")
save(prot_varpart, file=myfile, version=2)
