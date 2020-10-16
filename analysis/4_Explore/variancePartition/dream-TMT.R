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


#data(varPartDEdata)

# parallel processing params
params = BiocParallel::SnowParam(n_cores, "SOCK", progressbar=TRUE)
BiocParallel::register(params)

## munge input protein data, annotate with geno, subject and biofraction
geno <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
fraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
subject <- as.numeric(interaction(msstats_prot$Mixture,geno))

msstats_prot$Genotype <- geno
msstats_prot$Subject <- subject
msstats_prot$BioFraction <- fraction

## cast the data into a matrix
form = formula("Protein ~ Mixture + Channel + Genotype + BioFraction + Subject")
dm = msstats_prot %>% 
	reshape2::dcast(form, value.var="Abundance") %>% 
	as.data.table() %>% 
	na.omit() %>% 
	as.matrix(rownames="Protein") 

# collect sample metadata
metadata = as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(metadata) <- c("Mixture", "Channel", "Genotype", "BioFraction", "Subject")
rownames(metadata) <- colnames(dm)

# create dge object, perform TMM normalization, and subset for speed
dge = edgeR::DGEList(dm[c(1:10),])

# the model to be fit:
fx = formula("~ 0 + Genotype + BioFraction + (1|Subject)")

# do the repeated measures bit
dream_dge = variancePartition::voomWithDreamWeights(dge, fx, metadata)

# get contrast matrix for given coefficients
L = variancePartition::getContrast(dream_dge, fx, metadata, c("GenotypeMutant", "GenotypeControl"))

# fit DREAM model for each gene
dream_fit = variancePartition::dream(dream_dge, fx, metadata, L)

# extract results from first contrast

# get results
data(gene_map)
x = topTable(dream_fit, coef='GenotypeMutant', number=5 )
x = tibble::add_column(x,Protein=rownames(x),.before=1)
x = tibble::add_column(x,"Symbol" = gene_map$symbol[match(rownames(x),gene_map$uniprot)],.after="Protein")
