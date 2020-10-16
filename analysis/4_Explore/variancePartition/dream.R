#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root)

suppressPackageStartupMessages({
  library(dplyr)
  library(variancePartition)
})

data(varPartDEdata)

# filter genes based on threshold min = 0.1
idx = rowSums(cpm(countMatrix) > 0.1) >= 5

# create dge object, perform TMM normalization, and subset for speed
dge = edgeR::DGEList(countMatrix[idx,])
dge = edgeR::calcNormFactors(dge)
dge = dge[1:100,]

# parallel processing params
params = BiocParallel::SnowParam(4, "SOCK", progressbar=TRUE)
BiocParallel::register(params)

# the model to be fit:
fx = formula("~ Disease + (1|Individual)")

# do the repeated measures bit
dream_list = variancePartition::voomWithDreamWeights(dge, fx, metadata)

# fit DREAM model for each gene
fm_dream = variancePartition::dream(dream_list, fx, metadata)

# get results
topTable(fm_dream, coef='Disease1', number=5 ) %>% knitr::kable()

## advanced - single contrast

fx = formula("~ 0 + DiseaseSubtype + Sex + (1|Individual)")

# get contrast matrix for given coefficients
L = variancePartition::getContrast(dream_list, fx, metadata, c("DiseaseSubtype2", "DiseaseSubtype1"))

#plotContrasts(L) 
fit = variancePartition::dream(dream_list, fx, metadata, L)

colnames(fit)

# extract results from first contrast
topTable( fit, coef="L1", number=3 )


## define multiple contrasts
fx <- formula("~ 0 + DiseaseSubtype + Sex + (1|Individual)")
L1 = getContrast(dream_list, fx, metadata, c("DiseaseSubtype2", "DiseaseSubtype1"))
L2 = getContrast(dream_list, fx, metadata, c("DiseaseSubtype1", "DiseaseSubtype0"))
L = cbind(L1, L2)

#plotContrasts(L)

# fit the model
fit = dream(dream_list, fx, metadata, L)

# examine the results
topTable(fit, coef="L1", number=5) %>% knitr::kable()

# save data
#data.table::fwrite(countMatrix,"countMatrix.csv")
#data.table::fwrite(metadata,"metadata.csv")

L3 = c(1, -1/2, -1/2, 0)
Lall = cbind(L, data.frame(L3 = L3))

# joint hypothesis test!
topTable(fit, coef=c("DiseaseSubtype2", "DiseaseSubtype1"), number=3 )

# small dataset
fitmmKR = dream(dream_list, fx, metadata, ddf="Kenward-Roger")
