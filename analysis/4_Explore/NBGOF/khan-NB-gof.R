#!/usr/bin/env Rscript

# title: analysis of IOVS mouse lens data
# author: phillip wilmart
# description: analysis of TMT proteomics with edgeR

## input -------------------------------------------------------------
# * Kahn et al., Supplemental Table S01
data_url <- "https://bit.ly/2ZPu86T" # pwilmart/IRS_normalization/data.csv

# NOTE: data is TMT MS3

## options for NB GOF simulations:
nsim <- 999 # number of simulations to run
conf.env <- 0.95 # confidence interval
n_rows <- 500 # the number of genes to randomly sample and analyze
nthreads <- parallel::detectCores() - 1 # num of cores for parallel processing

# data source:
#	Khan, Shahid Y., et al. "Proteome Profiling of
#	Developing Murine Lens Through Mass Spectrometry."
# 	Investigative Ophthalmology & Visual Science 59.1
#	(2018): 100-107.

## Misc functions -------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# get the repository's root directory
	in_root <- function(h=here, dir=dpat) {
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE)))
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) {
		here <- dirname(here)
	}
	root <- here
	return(root)
}

# Prepare the R environment ---------------------------------------------------

# load renv
root <- getrd()
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(NBGOF)
	library(limma)
	library(edgeR)
	library(dplyr)
	library(data.table)
})

## load the data ----------------------------------------------------

#data_start <- read_csv("iovs-58-13-55_s01.csv")
raw_data <- data.table::fread(data_url)

# clean-up column names for numeric data
idy <- grep("Reporter",colnames(raw_data))
new_names <- sapply(strsplit(colnames(raw_data)[idy]," "),tail,1)
colnames(raw_data) <- c("Symbol","Accession",new_names)

# tidy up the data
df <- reshape2::melt(raw_data,
		     id.vars=c("Symbol","Accession"),
		     value.name="Intensity",variable.name="Sample") %>%
	as.data.table()

# Add Condition and Experiment annotations based on Sample name
df$Condition <- sapply(strsplit(as.character(df$Sample),"_"),"[",1)
df$Experiment <- sapply(strsplit(as.character(df$Sample),"_"),"[",2)

# collect a subset of samples cooresponding to one biological group (Condition)
all_groups <- unique(df$Condition)
subdf <- df %>% filter(Condition == sample(all_groups,1))

# cast the data into a matrix
dm <- subdf %>% dcast(Accession ~ Sample,value.var = "Intensity") %>%
	as.data.table() %>%
	as.matrix(rownames="Accession") %>%
	na.omit() # drop NA

# collect subset of random rows
idx <- sample(nrow(dm),n_rows)
subdm <- dm[idx,]

# create a DGEList object with our data
dge <- DGEList(counts = subdm, group =as.factor(rep(1,3)))

# TMM normalization
dge <- calcNormFactors(dge)

# estimate dispersion
dge <- estimateDisp(dge,design=as.matrix(rep(1,3)))

message(c("edgeR dispersion estimates: ",
	      paste(names(dge)[grep("dispersion",names(dge))],collapse=", ")))
#dge$common.dispersion
#dge$tagwise.dispersion -- I believe this is synonmyous with 'Genewise'
#dge$trended.dispersion

# GOF tests for different dispersion models:
# edgeR models: Common, Trended, Tagwise-Common,Tagwise-Trend
n_sim = 999
models <- c("Common"="CoxReid","Trended"="auto", "Genewise"="auto")
results <- list()
for (mod in names(models)){
	message(paste("Simulating", mod,
		      "dispersion model."))
	gof <- nb.gof.m(counts=dge$counts,
			x=as.matrix(rep(1,3)),
			sim=n_sim,
		        model=mod,
		        method=models[mod],
			ncores=nthreads)
	summary(gof,conv.env=0.95,data.node="SWIP TMT")
	results[[mod]] <- gof
}

# save the restults
myfile <- file.path(root,"rdata","khan_dispersion_simulations.RData")
saveRDS(results,file=myfile)
