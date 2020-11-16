#!/usr/bin/env Rscript 

# title: MSstatsTMT
# description: analysis of intrafraction comparisons with MSstats
# author: twab

## Input:
root <- "~/projects/SwipProteomics"


## Options
FDR_alpha = 0.05

## prepare the working environment ---------------------------------------------

# load renv
renv::load(root)

# load functions in root/R and make data in root/data accessible
# library(SwipProteomics)
devtools::load_all(root)

# load data in root/data
data(gene_map)
data(msstats_prot)
data(msstats_results)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  #library(MSstats) # twesleyb/MSstats
  #library(MSstatsTMT) # twesleyb/MSstats
})

# clean-up the data
msstats_prot$Run <- NULL
msstats_prot$TechRepMixture <- NULL
msstats_prot$Channel <- as.character(msstats_prot$Channel)
msstats_prot$BioReplicate <- as.character(msstats_prot$BioReplicate)
msstats_prot$Condition <- as.character(msstats_prot$Condition)
msstats_prot$Mixture <- as.character(msstats_prot$Mixture)
msstats_prot$Genotype <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 1)
msstats_prot$BioFraction <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 2)

proteins <- unique(as.character(msstats_prot$Protein))

data(gene_map)

idx <- match(msstats_prot$Protein,gene_map$uniprot)
msstats_prot <- msstats_prot %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")

idx <- match(msstats_results$Protein,gene_map$uniprot)
msstats_results <- msstats_results %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")

# cast the data into a matrix
dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

n_miss <- apply(dm,1,function(x) sum(is.na(x)))

# we can only impute up to 50% missing. Assume that the proteins with spuriously
keep <- n_miss < 0.5 * ncol(dm)
subdm <- dm[keep,]

knn_data <- impute::impute.knn(subdm)
knn_dm <- knn_data$data
knn_df <- reshape2::melt(knn_dm)
colnames(knn_df) <- c("Protein","Sample","Abundance")
knn_prot <- knn_df %>% mutate(Sample = as.character(Sample)) %>% 
	mutate(Mixture = sapply(strsplit(Sample,"_"),"[",1)) %>%
	mutate(Genotype = sapply(strsplit(Sample,"_"),"[",2)) %>%
	mutate(BioFraction = sapply(strsplit(Sample,"_"),"[",3))

# merge with msstats_prot
idy <- c("Protein","Abundance", "Mixture","Genotype","BioFraction")
foo = knn_prot %>% left_join(msstats_prot, by=idy)
any(is.na(foo$Abundance))

dm <- foo %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")
any(is.na(dm))


FDR_alpha <- 0.05
idx <- msstats_results$"adj.pvalue" < FDR_alpha
sigprots <- unique(as.character(msstats_results$Protein)[idx])




## format msstats_prot and resulst for downstream analysis -------------------

# combine statistical results
msstats_results <- rbind(results1,results2)

# clean-up the data
msstats_prot$Run <- NULL
msstats_prot$TechRepMixture <- NULL
msstats_prot$Channel <- as.character(msstats_prot$Channel)
msstats_prot$BioReplicate <- as.character(msstats_prot$BioReplicate)
msstats_prot$Condition <- as.character(msstats_prot$Condition)
msstats_prot$Mixture <- as.character(msstats_prot$Mixture)
msstats_prot$Genotype <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 1)
msstats_prot$BioFraction <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 2)

# annotate with gene Symbols and Entrez ids
idx <- match(msstats_prot$Protein,gene_map$uniprot)
msstats_prot <- msstats_prot %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")

idx <- match(msstats_results$Protein,gene_map$uniprot)
msstats_results <- msstats_results %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")

# cast the data into a matrix
dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

# number of missing vals
n_miss <- apply(dm,1,function(x) sum(is.na(x)))

# we can only impute up to 50% missing. 
# what is up with spurious missing vals?
keep <- n_miss < 0.5 * ncol(dm)
subdm <- dm[keep,]

knn_data <- impute::impute.knn(subdm)
knn_dm <- knn_data$data
knn_df <- reshape2::melt(knn_dm)
colnames(knn_df) <- c("Protein","Sample","Abundance")
knn_prot <- knn_df %>% mutate(Sample = as.character(Sample)) %>% 
	mutate(Mixture = sapply(strsplit(Sample,"_"),"[",1)) %>%
	mutate(Genotype = sapply(strsplit(Sample,"_"),"[",2)) %>%
	mutate(BioFraction = sapply(strsplit(Sample,"_"),"[",3))

# merge with msstats_prot
idy <- c("Protein","Abundance", "Mixture","Genotype","BioFraction")
msstats_prot <- knn_prot %>% left_join(msstats_prot, by=idy)

dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

stopifnot(!any(is.na(dm)))

# collect sig prots
idx <- msstats_results$"adj.pvalue" < FDR_alpha
sigprots <- unique(as.character(msstats_results$Protein)[idx])

