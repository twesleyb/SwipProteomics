#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw

## ---- input

# project's root directory
root = "~/projects/SwipProteomics"

# FDR threshold for differential abundance
FDR_alpha = 0.05


## ---- prepare the R workspace

# prepare the R workspace for the analysis

# load required packages and functions
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(data.table) # for working with tables
	library(doParallel) # for parallel processing
})

# load project specific functions and data
devtools::load_all(root, quiet=TRUE)


# ---- functions

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.


# function to register parallel backend
startCluster <- function(n_cores = parallel::detectCores() - 1, quiet=TRUE){
  if (!quiet) { message("Utilizing ", n_cores, " parallel processors.") }
  doParallel::registerDoParallel(cores = n_cores)
} #EOF

# stop parallel processing
stopCluster <- function() { doParallel::stopImplicitCluster() } 


## ---- load the data

data(swip_tmt)
data(swip_gene_map)


## ---- statistical analysis for overall Mutant-Controll contrast

# working with SL + IRS data, we normalize each protein to its Intensity sum
# across all BioFractions (sum normalization, see LOPIT approaches, e.g. Geledaki et al.), 
# this is refered to as relative_Intensity or rel_Intensity here.

# loop through all proteins, fit the model:
message("fit: log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)")

fx <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)

# example, fit the model to WASHC4 aka SWIP
swip <- gene_map$uniprot[grepl("Washc4", gene_map$symbol)]
fm <- lmerTest::lmer(fx, swip_tmt %>% subset(Protein == swip))

# NOTE: the fit is singular bc there is very little variance attributed to
# Mixture, for example, see variancePartition and its approach, implemented here
# by extracting the variance components with getVariance. 
#fm0 <- lme4::lmer(log2(Intensity) ~ (1|Condition) + (1|Mixture), data = swip_tmt %>% subset(Protein == swip))
#vp <- getVariance(fm0)
#pve <- vp/sum(vp)
#pve %>% t() %>% knitr::kable()

# create a contast with getContrast
LT <- getContrast(fm,"Mutant","Control")

# assess contrast with lmerTestContrast
lmerTestContrast(fm,LT) %>%
	mutate(Contrast='Mutant-Control') %>%
	mutate(Protein = swip) %>%
        unique() %>% knitr::kable()


## ----- dopar loop to fill all protein models

## Blargh... this used to work fine, but now seems to fail on any machine.

## loop to fit all proteins
proteins <- unique(swip_tmt$Protein)
# skip calculating derivs to speed things up
lmer_control <- lme4::lmerControl(check.conv.singular="ignore", calc.derivs=FALSE)

## NOTE: this is computationally intensive

# startCluster()
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

fit_list <- foreach(prot = proteins, .packages=c("dplyr")) %dopar% {
  fm <- lmerTest::lmer(fx, swip_tmt %>% subset(Protein == prot), control=lmer_control)
  return(fm)
}
names(fit_list) <- proteins

# stop parallel processing
stopCluster()


## ---- sequential loop if dopar fails

# if parallel loop fails, try sequential loop using %do%
#fit_list <- list()
#pbar <- txtProgressBar(max=length(proteins),style=3)
#pbar <- foreach(prot = proteins, .packages=c("dplyr")) %do% {
#  fit_list[[prot]] <- lmerTest::lmer(fx, swip_tmt %>% subset(Protein == prot), control=lmer_control)
#  setTxtProgressBar(pbar, value=match(prot,proteins))
#  return(pbar)
#}


## ---- test overall Mutant-Control comparison

# NOTE: I cannot explain the following error:
# In mclapply(argsList, FUN, mc.preschedule = preschedule, mc.set.seed = set.seed,  :

# why has parallel processing stopped working?

# loop to test overall contrast
#startCluster()
#results_list <- foreach(prot = proteins, .export = c("fit_list","LT"), .packages=c("dplyr")) %dopar% {
#	fm <- fit_list[[prot]]
#	result <- lmerTestContrast(fm, LT) %>%
#  	     mutate(Contrast='Mutant-Control') %>%
#	     mutate(Protein = prot) %>%
#	     unique()
#  return(result)
#} #EOL
#names(results_list) <- proteins
#stopCluster()


## ---- sequential loop to assess overall Mutant-Control comparison

# loop to test overall contrast
results_list <- list()
pbar <- txtProgressBar(max=length(proteins), style=3)
pbar <- foreach(prot = proteins) %do% {
	fm <- fit_list[[prot]]
	result <- lmerTestContrast(fm, LT) %>%
  	     mutate(Contrast='Mutant-Control') %>%
	     mutate(Protein = prot) %>%
	     unique()
  results_list[[prot]] <- result
  setTxtProgressBar(pbar, value=match(prot,proteins))
  return(pbar)
} #EOL


## ---- collect statistical results

df <- dplyr::bind_rows(results_list) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular) %>%
	arrange(Pvalue)


# check washc* prots results
washc_prots <- gene_map$uniprot[grepl("Washc*",gene_map$symbol)]
df %>% filter(Protein %in% washc_prots) %>%
	mutate(FDR=formatC(FDR)) %>%
	mutate(Pvalue=formatC(Pvalue)) %>%
	select(-Padjust,-Protein,-Entrez) %>%
	arrange(Symbol) %>% 
	knitr::kable() 

# results for overall WT v MUT comparison
mut_wt_results <- df


## ---- intra-BioFraction statistical analysis

# NOTE: does parallel processing still work here?

# NOTE: NOPE, I cannot explain the following error...
# In mclapply(argsList, FUN, mc.preschedule = preschedule, mc.set.seed =
# set.seed,  :
# scheduled cores 1, 2, 18 did not deliver results, all values of the jobs will
# be affected

# all mut and wt conditions
biofractions <- c("F4","F5","F6","F7","F8","F9","F10")
mut <- paste("ConditionMutant",biofractions,sep=".")
wt <- paste("ConditionControl",biofractions,sep=".")

# dopar loop
#startCluster()
results_list <- foreach(prot = proteins) %do% {
  # for a given protein fit, 
  # loop to evaluate B=7 'WT v Mut' intra-BioFraction contrasts 
  fm <- fit_list[[prot]]
  B <- 7
  res_list <- list()
  for (i in seq(B)) {
	LT <- getContrast(fm, mut[i],wt[i])
	res <- lmerTestContrast(fm,LT)
	res_list[[res$Contrast]] <- res
  }
  # collect results for all intra-BioFraction comparisons
  prot_results <- do.call(rbind, res_list) %>% mutate(Protein = prot)
  return(prot_results)
} #EOL
#stopCluster()

# parallel:
# startCluster()
# list <- foreach(i = seq(n)) %dopar% {
#   # do work
#   return(result)
# }
# stopCluster()

# sequential:
# NOTE: just exchange dopar for do!
# list <- foreach(i = seq(n)) %do% {
#   # do work
#   return(result)
# }
 
# so could easily be toggled with a switch and option like
# parallel = TRUE.

## BUT: doParallel::registerDoParallel... is having touble...

## Is it because we are creating par connection w/in a function??


## ---- collect and clean-up results

# collect results for all intra-BioFraction comparisons
df <- dplyr::bind_rows(results_list) %>%
	group_by(Contrast) %>%
	mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(Symbol = gene_map$symbol[match(Protein,gene_map$uniprot)]) %>%
	mutate(Entrez = gene_map$entrez[match(Protein,gene_map$uniprot)]) %>%
	dplyr::select(Protein, Symbol, Entrez, Contrast, log2FC, percentControl,
		      SE, Tstatistic, Pvalue, DF, S2, FDR, Padjust, isSingular)

# collect as named list
results <- df %>% group_by(Contrast) %>% group_split()
names(results) <- sapply(results,function(x) unique(x$Contrast))

# shorten names
namen <- names(results)
shorter <- gsub("ConditionMutant\\.F[0-9]{1,2}-|ConditionControl\\.","",namen)
names(results) <- shorter

# sort and combine with overall Mutant-Control results
all_results <- results[biofractions]
class(all_results) <- "list"
all_results[["Mutant-Control"]] <- mut_wt_results

# sort
all_results <- lapply(all_results, function(x) x %>% arrange(Pvalue))

# summary of sig results
sapply(all_results,function(x) sum(x$FDR<FDR_alpha)) %>%
	t() %>% knitr::kable()

# all statistical results!
swip_results <-  do.call(rbind, all_results)

# sig_prots
temp_df <- swip_results %>% filter(Contrast=='Mutant-Control')
sig_prots <- unique(temp_df$Protein[temp_df$FDR<FDR_alpha])


## ----  save key results

# save results as rda
myfile <- file.path(root,"data","swip_results.rda")
save(swip_results,file=myfile,version=2)
message("saved: ", myfile)

# write as excel
myfile <- file.path(root,"tables","SWIP-lmerTest-TMT-Results.xlsx")
write_excel(all_results, myfile)
message("saved: ", myfile)

# final normalized protein in tidy format as rda object
myfile <- file.path(datadir,"swip_tmt.rda")
save(swip_tmt,file=myfile,version=2)
message("saved: ", myfile)

# save sig_prots
myfile <- file.path(datadir,"swip_sig_prots.rda")
save(sig_prots, file=myfile,version=2)
message("saved: ", myfile)
