#!/usr/bin/env Rscript
#### Analysis in MSstatsTMT

root <- "~
## Load MSstatsTMT package
library(MSstatsTMT)

## Read Proteome Discoverer PSM report (MS3-based quantification)
raw <- read.delim("data/161117_SILAC_HeLa_UPS1_TMT10_5Mixtures_3TechRep_UPSdB_Multiconsensus_PD22_Intensity_03_with_FDR_control_PSMs.txt")
annotation <- read.csv("data/SpikeIn5mix_PD_annotation.csv")
head(raw)
head(annotation)


##############################
## Extra filtering for this dataset :
## 1: remove the proteins overlapped between spiked-in proteins and background proteins

# extract the list of spiked-in proteins
ups.prot <- as.character(unique(raw[grepl("ups", raw$Protein.Accessions), "Protein.Accessions"]))
# extract the list of background proteins
bg.prot <- as.character(unique(raw[!grepl("ups", raw$Protein.Accessions), "Protein.Accessions"]))
# overlapped proteins between spiked-in proteins and background proteins
inter <- ups.prot[which(gsub("ups", "", ups.prot) %in% bg.prot)]; inter

# generate the list of proteins to remove
protein.remove <- c(inter, gsub("ups", "", inter))

# remove the overlapped proteins 
raw <- raw[!raw$Protein.Accessions %in% protein.remove,]

##############################
## Make MSstatsTMT required format
##############################
required.input <- PDtoMSstatsTMTFormat(input = raw, 
                              annotation = annotation,
                              which.proteinid = "Protein.Accessions" ## same as default
                              )
head(required.input)

## count the number of proteins
length(unique(required.input$ProteinName))  #  4812

write.csv(required.input, file = "data/ControlMixture_PD_MS3_PSM_file.csv", row.names = FALSE)

##############################
## Protein summarization
## including protein-level summarization and normalization
##############################
processed.quant <- proteinSummarization(required.input,
                                     method = "msstats", # same as default
                                     global_norm=TRUE, # global peptide normalization between channels
                                     reference_norm=TRUE, # local protein normalization based on refernce channels
                                     MBimpute = TRUE,
                                     maxQuantileforCensored = NULL,
                                     remove_norm_channel = TRUE, # remove empty channels
                                     remove_empty_channel = TRUE # remove norm channels
                                     )

write.csv(processed.quant, file='data/ControlMixture_PD_MS3_proteinAbundance_byMSstatsTMT.csv', row.names = FALSE)

##############################
## Data visualization
##############################
# dataProcessPlotsTMT(data.psm = required.input,
#                     data.summarization = processed.quant,
#                     width = 15, height = 6,
#                     type = 'ProfilePlot',
#                     address = 'ControlMixture_PD_MS3_')

##############################
## Model-based group comparison + moderated t test + adjust p-value
##############################
test.MSstatsTMT <- groupComparisonTMT(data = processed.quant, 
                               contrast.matrix = "pairwise",
                               moderated = TRUE, # do moderated t test
                               adj.method = "BH") # multiple comparison adjustment

##############################
## save the result
##############################
write.csv(test.MSstatsTMT, file='data/ControlMixture_PD_MS3_testResult_byMSstatsTMT.csv', row.names = FALSE)
