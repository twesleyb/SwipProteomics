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
data(gene_map)
data(msstats_prot)
data(msstats_results)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
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


# map uniprot to gene symbols and entrez ids
proteins <- unique(as.character(msstats_prot$Protein))
idx <- match(msstats_prot$Protein,gene_map$uniprot)
msstats_prot <- msstats_prot %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")


colnames(msstats_prot)

any(is.na(msstats_prot$Abundance))

# why are there na?
msstats_prot[is.na(msstats_prot$Abundance),]

filt_prot <- msstats_prot %>% filter(!is.na(Abundance))

any(is.na(filt_prot$Abundance))

# cast the data into a matrix
dm <- filt_prot %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

head(dm)


dim(dm)


any_missing <- apply(dm,1,function(x) any(is.na(x)))

sum(any_missing)
# 567 proteins with missingness

n_missing <- apply(dm,1,function(x) sum(is.na(x)))

# types of missingness
table(n_missing)


# we impute missing values for proteins 
# we cannot impute proteins with more than 50% missingness
to_drop <- apply(dm,1,function(x) sum(is.na(x))> 0.5 * ncol(dm))

sum(to_drop)

# proteins for which we impute missing values
imp_prot <- names(which(apply(dm[!to_drop,],1,function(x) any(is.na(x)))))

length(imp_prot)

# 66 proteins which we cannot work with
knn_data <- impute::impute.knn(dm[!to_drop,])

knn_dm <- knn_data$data

stopifnot(sum(is.na(knn_dm))==0)

knn_df <- reshape2::melt(knn_dm)

colnames(knn_df) <- c("Protein","Sample","Abundance")

knn_prot <- knn_df %>% mutate(Sample = as.character(Sample)) %>% 
	mutate(Mixture = sapply(strsplit(Sample,"_"),"[",1)) %>%
	mutate(Genotype = sapply(strsplit(Sample,"_"),"[",2)) %>%
	mutate(BioFraction = sapply(strsplit(Sample,"_"),"[",3))

# merge with msstats_prot meta data
idy <- intersect(colnames(knn_prot),colnames(filt_prot))
swip_tmt <- knn_prot %>% left_join(filt_prot, by=idy)

# check - no na
dm <- swip_tmt %>% 
	reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, 
			value.var="Abundance") %>%
	as.data.table() %>% as.matrix(rownames="Protein")
stopifnot(!any(is.na(dm)))


## format msstats_prot and resulst for downstream analysis -------------------

# drop NA
filt_results <- msstats_results %>% 
  # keep or remove SingleMeasurePerCondition?
	filter(!is.na(adj.pvalue)) %>% filter(is.na(issue)) %>% select(-issue)


colnames(filt_results)[colnames(filt_results) == "Label"] <- "Contrast"
colnames(filt_results)[colnames(filt_results) == "pvalue"] <- "Pvalue"
colnames(filt_results)[colnames(filt_results) == "adj.pvalue"] <- "FDR"

# summary
filt_results %>% group_by(Contrast) %>% summarize(nSig = sum(Padjust<FDR_alpha)) %>% knitr::kable()

# Bonferroni padjust
filt_results <- filt_results %>% group_by(Contrast) %>% mutate(Padjust = p.adjust(Pvalue,method="bonferroni"))

sigprots <- unique(as.character(msstats_results$Protein)[idx])

# annotate with gene Symbols and Entrez ids
idx <- match(filt_results$Protein,gene_map$uniprot)
filt_results <- filt_results %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")


colnames(filt_results)

results_list <- filt_results %>% group_by(Contrast) %>% group_split()

names(results_list) <- sapply(results_list,function(x) unique(x$Contrast))
namen <- names(results_list)
new_names <- gsub("Mutant.F[0-9]{1,2}-Control\\.","",namen)
names(results_list) <- new_names

# sort results list
idx <- c("F4","F5","F6","F7","F8","F9","F10","Mutant-Control")
results_list <- results_list[idx]

class(results_list) <- "list"

# save as excel document
myfile <- file.path(root,"tables","results.xlsx")
write_excel(results_list,myfile)
