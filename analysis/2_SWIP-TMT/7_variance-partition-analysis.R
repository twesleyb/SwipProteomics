#!/usr/bin/env Rscript

# prepare the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(swip)
data(partition)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(variancePartition)
})


## prepare the data -----------------------------------------------------------

# cast the data into a matrix.
fx <- formula(Protein ~ Mixture + Channel + BioFraction + Genotype + Subject)
dm <- msstats_prot %>%
	reshape2::dcast(fx, value.var= "Abundance") %>% 
	as.data.table() %>% na.omit() %>% as.data.table() %>% as.matrix(rownames="Protein")

# munge to create sample info from the dcast fx
info <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(info) <- strsplit(as.character(fx)[3]," \\+ ")[[1]]

# Create a formula specifying all covariates as mixed effects and do
# variancePartition -- calculate the percent variance explained.
form = formula(~ (1|Mixture) + (1|Channel) + (1|BioFraction) + (1|Genotype))

## analyze variance explained by covariates:
prot_varpart <- fitExtractVarPartModel(dm, form, info) %>% 
	as.data.frame() %>% as.data.table(keep.rownames="Protein")


## save ------------------------------------------------------------------------

myfile <- file.path(root,"data","prot_varpart.rda")
save(prot_varpart, file=myfile, version=2)

myfile <- file.path(root,"tables","S4_variancePartition_Results.xlsx")
write_excel(list("Variance Partition" = prot_varpart),myfile)
