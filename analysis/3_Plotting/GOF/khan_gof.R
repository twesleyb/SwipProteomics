#!/usr/bin/env Rscript

# title: analysis of IOVS mouse lens data
# author: phillip wilmart
# description: analysis of TMT proteomics with edgeR

## input -------------------------------------------------------------
# * Kahn et al., Supplemental Table S01
data_url <- "https://bit.ly/2ZPu86T" # pwilmart/IRS_normalization/data.csv

# data source:
#	Khan, Shahid Y., et al. "Proteome Profiling of
#	Developing Murine Lens Through Mass Spectrometry."
# 	Investigative Ophthalmology & Visual Science 59.1
#	(2018): 100-107.

## main -------------------------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(limma)
	library(edgeR)
	library(dplyr)
	library(data.table)
})

# download the data
raw_data <- data.table::fread(data_url)

# clean-up column names for numeric data
idy <- grep("Reporter",colnames(raw_data))
new_names <- sapply(strsplit(colnames(raw_data)[idy]," "),tail,1)
colnames(raw_data) <- c("Symbol","Accession",new_names)

# filter out proteins not seen in all three runs
data_no_na <- na.omit(raw_data)

# save the annotation columns (gene symbol and protein accession) for later and remove from data frame
annotate_df <- data_no_na[,1:2]
data_raw <- data_no_na[,3:20]
row.names(data_raw) <- annotate_df$Accession

# separate the TMT data by experiment
exp1_raw <- data_raw[,c(1:6)]
exp2_raw <- data_raw[,c(7:12)]
exp3_raw <- data_raw[,c(13:18)]

# figure out the global scaling value
target <- mean(c(colSums(exp1_raw),
		 colSums(exp2_raw),
		 colSums(exp3_raw)))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
# seems like a loop could be used here somehow...
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

# make working frame with row sums from each frame
irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl), rowSums(exp3_sl))
colnames(irs) <- c("sum1", "sum2", "sum3")

# get the geometric average intensity for each protein
irs$average <- apply(irs, 1, function(x) exp(mean(log(x))))

# compute the scaling factor vectors
irs$fac1 <- irs$average / irs$sum1
irs$fac2 <- irs$average / irs$sum2
irs$fac3 <- irs$average / irs$sum3

# make new data frame with normalized data
data_irs <- exp1_sl * irs$fac1
data_irs <- cbind(data_irs, exp2_sl * irs$fac2)
data_irs <- cbind(data_irs, exp3_sl * irs$fac3)

# set up the sample mapping
group <- factor(sapply(strsplit(colnames(data_irs),"_"),"[",1),
		levels = c("E15", "E18", "P0", "P3", "P6", "P9"))

# create a DGEList object with our data
y_sl <- DGEList(counts = data_irs, group = group)

y_sl$samples$batch <- factor(sapply(strsplit(rownames(y_sl$samples),"_"),"[",1),
		       levels = c("Set1","Set2","Set3"))

# perform TMM normalization
y_sl <- calcNormFactors(y_sl)

# estimate all dispersion metrics
y_sl <- estimateDisp(y_sl) 

## plot gof for three dispersion models

# plot common gof
myfile <- file.path(root,"figs","GOF","khan-common-gof.pdf")
pdf(file=myfile)
fit <- glmFit(y_sl,dispersion=y_sl$common.dispersion)
nbglm_gof <- gof(fit, plot=TRUE)
invisible(dev.off())
message("Outliers: ", sum(nbglm_gof$outlier))

# plot trended gof
myfile <- file.path(root,"figs","GOF","khan-trended-gof.pdf")
pdf(file=myfile)
fit <- glmFit(y_sl,dispersion=y_sl$trended.dispersion)
nbglm_gof <- gof(fit, plot=TRUE)
invisible(dev.off())
message("Outliers: ", sum(nbglm_gof$outlier))

# plot tagwise gof
myfile <- file.path(root,"figs","GOF","khan-tagwise-gof.pdf")
pdf(file=myfile)
fit <- glmFit(y_sl,dispersion=y_sl$tagwise.dispersion)
nbglm_gof <- gof(fit, plot=TRUE)
invisible(dev.off())
message("Outliers: ", sum(nbglm_gof$outlier))
