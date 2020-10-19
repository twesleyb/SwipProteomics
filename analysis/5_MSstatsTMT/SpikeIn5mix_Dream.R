#!/usr/bin/env Rscript

## input
n_cores = 23

## prepare the env
root <- "~/projects/SwipProteomics"; renv::load(root)
devtools::load_all(root)

# load pd_annotation
myfile <- file.path(root,"rdata","ms3_pd_annotation.rda")
load(myfile)

# load msstats_results
myfile <- file.path(root,"rdata","ms3_msstats_results.rda")
load(myfile)

# load msstats_prot
myfile <- file.path(root,"rdata","ms3_msstats_prot.rda")
load(myfile)

# parallel processing
BioCParallel(register(BiocParallel::SnowParam(n_cores)))

## imports
library(dplyr)
library(MSstatsTMT)
library(variancePartition)

# cast the data into a matrix, drop NA
dm <- msstats_prot %>% 
	reshape2::dcast(Protein ~ Mixture + TechRepMixture + Channel + Condition + BioReplicate, value.var = "Abundance") %>%
	as.data.table() %>% na.omit() %>% as.matrix(rownames="Protein")

# create sample meta info
samples <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(samples) <- c("Mixture", "TechRepMixture", "Channel", 
		       "Condition", "BioReplicate")
rownames(samples) <- colnames(dm)

# Abundance ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture) + Condition
form <- ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture) + Condition

# fit the model
fit <- dream(dm, form, samples)

# test a contrast
L1 <- getContrast(dm, form, samples, c("Condition0.5","Condition1"))
fit <- dream(dm, form, samples, L1)

# collect results
results_df <- topTable(fit, coef="L1",number=Inf)

results_df %>% head() %>% knitr::kable()

# msstats_results for the given contrast:
results_list <- msstats_results %>% group_by(Label) %>% group_split()
names(results_list) <- sapply(results_list,function(x) unique(x$Label))

df = results_list[["0.5-1"]] %>% arrange(pvalue) 
df %>% head() %>% knitr::kable()

#names(results_list)
