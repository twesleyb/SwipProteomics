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
form <- formula("~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Subject)")

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

# save
myfile <- file.path(root,"rdata","dream_results.csv")
fwrite(results_df,file=myfile)

# status
message("\nNumber of differentially abundant proteins: ", 
	sum(results_df$adj.P.Val < 0.05))
results_df %>% filter(adj.P.Val < 0.05) %>% knitr::kable(row.names=FALSE)

################################################################
## what about BioFraction versus all else?

# function
annotate <- function(x) {
	protein <- rownames(x)
	x = tibble::add_column(x,protein,.before=1)
	idx <- match(protein, gene_map$uniprot)
	symbol = gene_map$symbol[idx]
	x = tibble::add_column(x,symbol,.after=1)
	return(x)
}


## the model to be fit:
form <- formula("~ 0 + BioFraction + Genotype + (1|Mixture) + (1|Subject)")

# munge up some contrasts
contrasts <- c("F4","F5","F6","F7","F8","F9","F10")
L0 = c(BioFractionF4 = NA, BioFractionF5 = NA, # F4 
       BioFractionF6 = NA, BioFractionF7 = NA, 
       BioFractionF8 = NA, BioFractionF9 = NA, 
       BioFractionF10 = NA, GenotypeMutant = 0)
contrasts_list <- list()
for (i in c(1:length(contrasts))) {
  L = L0
  L[which(grepl(contrasts[i],names(L)))] <- 1
  L[is.na(L)] <- -1/6
  contrasts_list[[i]] <- L
}
names(contrasts_list) <- contrasts

L <- do.call(cbind, contrasts_list)

# fit the model
fit <- dream(dm, form, samples, L)

# collect the results
results_list <- lapply(colnames(L),function(x) topTable(fit,coef=x,number="Inf"))
names(results_list) <- contrasts

# annotate with gene ids
results_list <- lapply(results_list, annotate)


write_excel(results_list,"dream-BioFraction-vs-else.xlsx")
