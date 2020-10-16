#!/usr/bin/env Rscript

## options
n_cores = 23

## input
root = "~/projects/SwipProteomics"

## prepare the env

renv::load(root)
devtools::load_all()

# local data
data(msstats_prot)

# global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(data.table)
  library(variancePartition)
})

# parallel processing params
params = BiocParallel::SnowParam(n_cores, "SOCK", progressbar=TRUE)
BiocParallel::register(params)

## munge input protein data, annotate with geno, subject and biofraction
geno <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
subject <- as.numeric(interaction(msstats_prot$Mixture,geno))

msstats_prot$Genotype <- geno
msstats_prot$Subject <- subject
msstats_prot$BioFraction <- biofraction

## cast the data into a matrix
form = formula("Protein ~ Mixture + Channel + Genotype + BioFraction + Subject")
dm = msstats_prot %>% 
	reshape2::dcast(form, value.var="Abundance") %>% 
	as.data.table() %>% 
	na.omit() %>% 
	as.matrix(rownames="Protein") 

# collect sample metadata
metadata = as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(metadata) <- c("Mixture", "Channel", "Genotype", 
			"BioFraction", "Subject")
rownames(metadata) <- colnames(dm)

# create dge object, perform TMM normalization, and subset for speed
data(swip)

#dge = edgeR::DGEList(dm[c(1:10),])
dge = edgeR::DGEList(dm[c(swip,sample(rownames(dm),9)),])

# the model to be fit:
fx = formula("~ Genotype + (1|Subject) + BioFraction")

# do the repeated measures bit
dream_dge = variancePartition::voomWithDreamWeights(dge, fx, metadata)

# get contrast matrix for given coefficients
#L = variancePartition::getContrast(dream_dge, fx, metadata, "GenotypeControl")

# fit DREAM model for each gene
# NOTE: you don't have to specify a contrast L
# This actually might be the correct contrast!
dream_fit = variancePartition::dream(dream_dge, fx, metadata)

# get results
data(gene_map)
x = topTable(dream_fit, coef='GenotypeMutant', number=5 )
x = tibble::add_column(x,Protein=rownames(x),.before=1)
idx = match(rownames(x),gene_map$uniprot)
x = tibble::add_column(x,"Symbol" = gene_map$symbol[idx], .after="Protein")
knitr::kable(x)
