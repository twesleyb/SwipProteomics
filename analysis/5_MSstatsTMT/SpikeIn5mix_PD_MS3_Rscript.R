#!/usr/bin/env Rscript

# title: analysis using MSstatsTMT

root <- "~/projects/SwipProteomics"
renv::load(root)

## Load MSstatsTMT package
library(MSstatsTMT)

## Read Proteome Discoverer PSM report (MS3-based quantification)
myfile <- file.path(root,"rdata","SpikeIn5mix_MS3_PD_PSM.txt")
raw <- read.delim(myfile)

## load annotations
myfile <- file.path(root,"data","SpikeIn5mix_PD_annotation.csv")
annotation <- read.csv(myfile)


## Extra filtering for this dataset :
## 1: remove the proteins overlapped between spiked-in proteins and background
## proteins

# extract the list of spiked-in proteins
ups.prot <- as.character(unique(raw[grepl("ups", raw$Protein.Accessions), 
				"Protein.Accessions"]))
# extract the list of background proteins
bg.prot <- as.character(unique(raw[!grepl("ups", raw$Protein.Accessions), 
			       "Protein.Accessions"]))
# overlapped proteins between spiked-in proteins and background proteins
inter <- ups.prot[which(gsub("ups", "", ups.prot) %in% bg.prot)]; inter

# generate the list of proteins to remove
protein.remove <- c(inter, gsub("ups", "", inter))

# remove the overlapped proteins 
raw <- raw[!raw$Protein.Accessions %in% protein.remove,]


## Make MSstatsTMT required format
required.input <- MSstatsTMT::PDtoMSstatsTMTFormat(input = raw, annotation)


## Protein summarization
## including protein-level summarization and normalization
processed.quant <- MSstatsTMT::proteinSummarization(required.input,
                                     method = "msstats", 
                                     global_norm=TRUE, 
                                     reference_norm=TRUE, 
                                     MBimpute = TRUE,
                                     maxQuantileforCensored = NULL,
                                     remove_norm_channel = TRUE, 
                                     remove_empty_channel = TRUE,
				     clusters=23
                                     )


## Model-based group comparison + moderated t test + adjust p-value
test.MSstatsTMT <- MSstatsTMT::groupComparisonTMT(data = processed.quant, 
                               contrast.matrix = "pairwise",
                               moderated = TRUE, # do moderated t test
                               adj.method = "BH") 


## save the result
myfile <- file.path(root,"rdata",
		    "ControlMixture_PD_MS3_testResult_byMSstatsTMT.csv")
write.csv(test.MSstatsTMT, file=myfile, row.names=FALSE)
