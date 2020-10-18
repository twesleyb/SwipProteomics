#!/usr/bin/env Rscript

## options
n_cores = 23

## input
root = "~/projects/SwipProteomics"

## prepare the env
renv::load(root)
devtools::load_all()

# local data
data(swip)
data(gene_map)
data(msstats_prot)
data(leidenalg_partition)

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
metadata <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(metadata) <- c("Mixture", "Channel", "Genotype", 
			"BioFraction", "Subject")
rownames(metadata) <- colnames(dm)

knitr::kable(metadata)


# the model to be fit:
form <- formula("~ BioFraction + (1|Subject) + Genotype")

## create dge object
dge <- edgeR::DGEList(dm)

# do the repeated measures bit
dream_dge <- variancePartition::voomWithDreamWeights(dge, form, metadata)

# names(dream_dge)
# [1] "targets" "E"       "weights" 

# create a contrast to be tested:
# NOTE: simple contrasts are tested by default
#L <- getContrast(dream_dge, form, metadata, "GenotypeMutant")

# the variance for each protein parititioned
#fx = ~ (1|Mixture) + (1|Channel) + (1|BioFraction) + (1|Subject) + (1|Genotype)
#vp = fitExtractVarPartModel(dream_dge, fx, metadata)
#head(vp) # typically, BioFraction explains a vast majority of the variance
#vp[swip,] # for some proteins Genotype explains the majority of variance
#plotVarPart( sortCols(vp))

# fit DREAM model for each protein with comparisons specified by contrasts L
#dream_fit <- variancePartition::dream(dream_dge, form, metadata, L)
dream_fit <- variancePartition::dream(dream_dge, form, metadata)

# get results
results <- topTable(dream_fit, number=Inf)

# clean-up results
results <- tibble::add_column(results, Protein=rownames(results),.before=1)
idx <- match(results$Protein,gene_map$uniprot)
symbols <- gene_map$symbol[idx]
results <- tibble::add_column(results,Symbol=symbols,.after="Protein")

# inspect top proteins
results %>% arrange(P.Value) %>% knitr::kable()
