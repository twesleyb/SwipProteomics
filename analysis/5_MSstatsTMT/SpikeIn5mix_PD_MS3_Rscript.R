#!/usr/bin/env Rscript

## prepare the env
root <- "~/projects/SwipProteomics"; renv::load(root)

## imports
library(dplyr)
library(MSstatsTMT)

## Read Proteome Discoverer PSM report (MS3-based quantification)
myfile <- file.path(root,"rdata","SpikeIn5mix_MS3_PD_PSM.txt")
raw <- read.delim(myfile)

## load sample annotations
myfile <- file.path(root,"data","SpikeIn5mix_PD_annotation.csv")
annotation <- read.csv(myfile)


## Extra filtering 
# remove the proteins overlapped between spiked-in proteins and background
# proteins

# extract the list of spiked-in proteins
ups.prot <- as.character(unique(raw[grepl("ups", raw$Protein.Accessions), 
				"Protein.Accessions"]))

# extract the list of background proteins
bg.prot <- as.character(unique(raw[!grepl("ups", raw$Protein.Accessions), 
			       "Protein.Accessions"]))

# overlapped proteins between spiked-in proteins and background proteins
inter <- ups.prot[which(gsub("ups", "", ups.prot) %in% bg.prot)]

# generate the list of proteins to remove
protein.remove <- c(inter, gsub("ups", "", inter))

# remove the overlapped proteins 
raw <- raw[!raw$Protein.Accessions %in% protein.remove,]


## Make MSstatsTMT required format
required.input <- suppressMessages({
  MSstatsTMT::PDtoMSstatsTMTFormat(input = raw, annotation)
})


## Protein summarization
## including protein-level summarization and normalization
msstats_prot <- suppressMessages({
	MSstatsTMT::proteinSummarization(required.input,
                                     method = "msstats", 
                                     global_norm=TRUE, 
                                     reference_norm=TRUE, 
                                     MBimpute = TRUE,
                                     maxQuantileforCensored = NULL,
                                     remove_norm_channel = TRUE, 
                                     remove_empty_channel = TRUE,
				     clusters=23)
})


## Model-based group comparison + moderated t test + adjust p-value
msstats_results <- suppressMessages({
	MSstatsTMT::groupComparisonTMT(data = msstats_prot, 
                               contrast.matrix = "pairwise",
                               moderated = TRUE, # do moderated t test
                               adj.method = "BH") 
})


# check which model was fit to the data:
protein <- sample(unique(as.character(msstats_prot$Protein)),1)

# Abundance ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture) + Condition
fit <- MSstatsTMT::checkDesign(msstats_prot,protein)

summary(fit)

#r.squaredGLMM.merMod(fit)
#           R2m       R2c 
#[1,] 0.02981312 0.3596938

# examine sig
msstats_results %>% 
	group_by(Label) %>% 
	summarize(nsig=sum(adj.pvalue < 0.05, na.rm=TRUE)) %>% 
	knitr::kable()


## save the results -----------------------------------------------------------

pd_annotation <- annotation
myfile <- file.path(root,"rdata","ms3_pd_annotation.rda")
save(pd_annotation, file=myfile,version=2)

myfile <- file.path(root,"rdata","ms3_msstats_results.rda")
save(msstats_results,file=myfile,version=2)

myfile <- file.path(root,"rdata","ms3_msstats_prot.rda")
save(msstats_prot,file=myfile,version=2)
