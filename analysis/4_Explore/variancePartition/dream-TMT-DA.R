#!/usr/bin/env Rscript

## options
n_cores = 23

## input
root = "~/projects/SwipProteomics"

## prepare the env
renv::load(root)
devtools::load_all()

# local data
data(gene_map) 
data(msstats_prot)

# global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(data.table)
  library(variancePartition)
})


## parallel processing params
params <- BiocParallel::SnowParam(n_cores, "SOCK", progressbar=TRUE)
BiocParallel::register(params)


## munge input protein data, annotate with geno, subject and biofraction
geno <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
msstats_prot$Genotype <- factor(geno,levels=c("Control","Mutant"))
msstats_prot$BioFraction <- factor(biofraction, 
				   levels=c("F4","F5","F6","F7","F8","F9","F10"))
subject <- as.numeric(interaction(msstats_prot$Mixture,msstats_prot$Genotype))
msstats_prot$Subject <- subject


## cast the data into a matrix
fx <- formula("Protein ~ Mixture + Channel + Genotype + BioFraction + Subject")
dm <- msstats_prot %>% 
	reshape2::dcast(fx, value.var="Abundance") %>% 
	as.data.table() %>% 
	na.omit() %>% 
	as.matrix(rownames="Protein") 

# subset
#dm <- dm[sample(nrow(dm),10),]

## collect sample metadata
samples <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(samples) <- c("Mixture", "Channel", "Genotype", 
			"BioFraction", "Subject")
rownames(samples) <- colnames(dm)
samples$Genotype <- factor(samples$Genotype)
samples$Subject <- factor(samples$Subject)
samples$BioFraction <- factor(samples$BioFraction,
			      levels=c("F4","F5","F6","F7","F8","F9", "F10"))

levels(samples$Genotype)
levels(samples$Subject)
# There are six mice (subjects)--3x Control and 3x Mutant.

levels(samples$BioFraction)
# We measured 7 subcellular fractions (BioFraction) per mouse.

# the model to be fit:
# NOTE: Q1 to include an intercept (~ 1 + ...) or not?
form <- formula("~ 1 + Genotype + BioFraction + (1|Subject)")
#form <- formula("~ Genotype + BioFraction + (1|Subject)")
#form <- formula("~ 0 + Genotype + BioFraction + (1|Subject)")


## create dge object
dge <- edgeR::DGEList(dm)


## do the repeated measures bit with variancePartition (voom)
dream_dge <- variancePartition::voomWithDreamWeights(dge, form, samples)

# names(dream_dge)
# [1] "targets" "E"       "weights" 


## create a contrast to be tested:
# NOTE: simple contrasts are tested by default
L <- getContrast(dream_dge, form, samples, c("GenotypeMutant","(Intercept)"))

# check the design:
plotContrasts(L)

## fit DREAM model for each protein with comparisons specified by contrasts L
dream_fit <- variancePartition::dream(dream_dge, form, samples, L)

# Warning message:
# In variancePartition::dream(dream_dge, form, samples, L) :
# Contrasts with only a single non-zero term are already evaluated by default.


## get results
results <- topTable(dream_fit, coef="L1", number=Inf)

# clean-up results
results <- tibble::add_column(results, Protein=rownames(results),.before=1)
idx <- match(results$Protein,gene_map$uniprot)
symbols <- gene_map$symbol[idx]
results <- tibble::add_column(results,Symbol=symbols,.after="Protein")

# inspect top proteins
results %>% arrange(P.Value) %>% head() %>% knitr::kable()
