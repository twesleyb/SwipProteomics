#!/usr/bin/env Rscript

## options
n_cores = 23


## prepare the env
root = "~/projects/SwipProteomics"; renv::load(root)
suppressWarnings({ devtools::load_all() })


## parallel processing params
params <- BiocParallel::SnowParam(n_cores, "SOCK", progressbar=TRUE)
BiocParallel::register(params)


## local data
data(swip)
data(gene_map) 
data(msstats_prot)


## global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(variancePartition)
})


## cast the data into a matrix
fx <- formula("Protein ~ Mixture + Channel + Genotype + BioFraction + Subject")
dm <- msstats_prot %>% 
	reshape2::dcast(fx, value.var="Abundance") %>% 
	as.data.table() %>% 
	na.omit() %>% 
	as.matrix(rownames="Protein") 


## collect sample metadata
samples <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(samples) <- c("Mixture", "Channel", "Genotype", 
			"BioFraction", "Subject")

# important! rownames of sample metadata must be colnames of expression data
rownames(samples) <- colnames(dm)
samples$Genotype <- factor(samples$Genotype)
samples$Subject <- factor(samples$Subject)
samples$BioFraction <- factor(samples$BioFraction,
			      levels=c("F4","F5","F6","F7","F8","F9", "F10"))

# There are six mice (subjects)--3x Control and 3x Mutant.
# We measured 7 subcellular fractions (BioFraction) per mouse.


## the model to be fit:
form <- formula("~ 0 + Genotype + BioFraction + (1|Subject)")

#L1 = getContrast(subdm, form, samples, c("GenotypeMutant","GenotypeControl"))
L1 = c(GenotypeControl = -1, GenotypeMutant = 1, 
       BioFractionF5 = 0, BioFractionF6 = 0, BioFractionF7 = 0, 
       BioFractionF8 = 0, BioFractionF9 = 0, BioFractionF10 = 0)

fit <- dream(dm, form, samples, L1)


## collect results
results_df <- topTable(fit, coef="L1",number=Inf)

Protein <- rownames(results_df)
Symbol <- gene_map$symbol[match(Protein,gene_map$uniprot)]
results_df <- tibble::add_column(results_df,Protein,.before=1)
results_df <- tibble::add_column(results_df,Symbol,.after=1)

knitr::kable(results_df, row.names=FALSE)
