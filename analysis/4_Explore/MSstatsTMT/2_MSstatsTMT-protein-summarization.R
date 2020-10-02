#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

# FIXME: catch/suppress warnings about singular fit?
# FIXME: exchange pbar for messages

## Input ----------------------------------------------------------------------
# specify the projects root directory:
root = "~/projects/SwipProteomics"

# the PSM data is in root/data:
input_dir = "data/PSM.zip"

# PSM.zip contains:
input_psm = "PSM/5359_PSM_Report.xlsx"
input_samples = "PSM/5359_Sample_Report.xlsx"

## Options --------------------------------------------------------------------
load_rda = FALSE # load saved objects to speed things up
save_rda = FALSE # save key R objects?
n_prots = 10 # subset of proteins to analyze for speed

## Output ---------------------------------------------------------------------
# * several intermediate datasets in root/rdata
# * the MSstatsTMT statistical results for intrafraction comparisons saved as
#   an excel workbook in root/tables

## DESCRIPTION:

# * fit appropriate model
# * get some basic stats from fit with anova(), coeff()
# do this thing rhoa to get fixEffs, sigma thopt
# calculate asymptotic var-covar matrix for theata and sigma
#rho$fixEffs <- lme4::fixef(rho$model)
#rho$sigma <- stats::sigma(rho$model)
#rho$thopt <- lme4::getME(rho$model,"theta")

# given rho, calcApvar

# calcApvar <- function(rho){
# dd <- .devfunTheta(rho$model) # see [10]
# #update(fm, devFunOnly = TRUE)
# 
#   # calculate 2x2 hessian matrix
#   h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma)) # see [11]
#   #       [,1]        [,2]
#   # [1,]  18.03627   148.6282
#   # [2,] 148.62825 17146.7792
#   # calculate Choleski factorization of h (decomposition)
#   ch <- try(base::chol(h), silent = TRUE)
#   #       [,1]        [,2]
#   # [1,] 4.246914  34.99677
#   # [2,] 0.000000 126.18243
#   if (inherits(ch, "try-error")) {
#     return(rho)
#   }
#   # calculate the inverse of the Choleski matrix
#   A <- 2 * base::chol2inv(ch)
#   #       [,1]        [,2]
#   # [1,]  0.119417491 -0.0010351106
#   # [2,] -0.001035111  0.0001256123
#   # calculate eigenvalue of hessian matrix h (spectral decomposition)
#   eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
#   # [1] 17148.06878    16.74671 
#   # tol ~ sqrt(.Machine$double.eps)
#   isposA <- !(min(eigval) < sqrt(.Machine$double.eps))
#   if (!isposA) {
#     warning("Asymptotic covariance matrix A is not positive!")
#   }
#   # return inverse of cholesky matrix of hessian d
#   # NOTE: important bits:
#   # looks like .myhess is some sort of wrapper around hessian, using devFunTheta
#   #  >>>  dd <- .devfunTheta(rho$model) 
#   #  >>>  h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma))
#   #  >>>  ch <- try(base::chol(h), silent = TRUE)
#   #  >>>  A <- 2 * base::chol2inv(ch)
#   return(A)
# }


## FUNCTIONS -----------------------------------------------------------------

squote <- function(string) { 
	# wrap a string in singular quotation marks (')
	return(paste0("'",string,"'"))
}


mkdir <- function(...,warn=TRUE,report=FALSE) {
	# create a new directory
	newdir <- file.path(...)
	if (warn & dir.exists(newdir)) {
		warning("dir exists")
	} else if (!dir.exists(newdir)) { 
		dir.create(newdir)
		if (report) {
		message(paste("created",newdir))
		}
	}
}


fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}


munge1 <- function(x) { # x = samples$ConditionFraction
	# function to reformat the data
	# extract sample 'Condition' annotation from ConditionFraction
	paste0("F",as.numeric(sapply(strsplit(x,"Control|Mutant|SPQC"),"[",2)))
}


munge2 <- function(x) { # x = samples$ConditionFraction 
	# function to reformat the data
	# extract sample 'Fraction' annotation from ConditionFraction
	sapply(strsplit(x,"[0-9]{1,2}"),"[",1)
}


munge3 <- function(x) { # x = samples$Experiment
	# coerce Experiment to replicate ID = R#
	paste0("R",as.numeric(as.factor(samples$Experiment))) 
}


reformat_cols <- function(raw_pd) {
	# make columns look like MSstats by replacing special characters with '.'
	# replace special characters in column names with "."
	chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
	new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
	colnames(raw_pd) <- new_names
	# add 'X' if starts column name starts with ".."
	colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))
	# return the reformatted data
	return(raw_pd)
}


## Prepare the working environment --------------------------------------------

# project directories
datadir <- file.path(root,"data"); mkdir(datadir,warn=FALSE)
rdatdir <- file.path(root,"rdata"); mkdir(rdatdir,warn=FALSE)
downdir <- file.path(root,"downloads"); mkdir(downdir,warn=FALSE)


# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	suppressWarnings({ library(getPPIs) }) # FIXME: annoying warnings!
	library(data.table)
	library(MSstatsTMT) # twesleyb/MSstatsTMT
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## load preprocessed data ------------------------------------------------------

load(file.path(root,"rdata","data_pd.rda"))

# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
if (load_rda) {
  load(file.path(root,"rdata","msstats_prot.rda"))
  message("Loaded 'msstats_prot'")
} else {
  message("\nPerforming normalization and protein sumamrization.")
  msstats_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)
}


if (save_rda) {
  # save to file
  myfile <- file.path(rdatdir, "msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message(paste("\nSaved",myfile))
}


## Remove protein outliers? -----------------------------------------------------
	
# Subset for speed:
rand_prots <- sample(unique(msstats_prot$Protein),n_prots)

filt_prot <- msstats_prot %>% filter(Protein %in% rand_prots)


## Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
# FIXME: parallelize!

foo = capture.output({
	all_results <- groupComparisonTMT(filt_prot)	
})

# there are 9x calls per protein to thpars ... 
# seems to be a call to update(fm) 
# maybe from the fit of P42225:
#[1] "numeric"
#Run.(Intercept)           sigma 
#     0.54006363      0.08232073 

#load(file.path(rdatdir,"all_results.rda"))

# save as rda
if (save_rda) {
	myfile <- file.path(root,"rdata","all_results.rda")
	save(all_results,file=myfile,version=2)
	message(paste("\nsaved",squote(basename(myfile)),"in",
		      squote(dirname(myfile))))
}

quit()

## clean-up MSstats results ----------------------------------------------------

# collect relevant contrasts
control_groups <- paste("Control",paste0("F",seq(4,10)),sep=".")
mutant_groups <- paste("Mutant",paste0("F",seq(4,10)),sep=".")
contrasts <- paste(control_groups,mutant_groups,sep="-")
filt_results <- all_results %>% filter(Label %in% contrasts)

# Symbol and Entrez to the data
idx <- match(filt_results$Protein,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
filt_results <- tibble::add_column(filt_results,Entrez,.after="Protein")
filt_results <- tibble::add_column(filt_results,Symbol,.after="Entrez")

# Can we add protein level data to  stats and save as tidy object?
# format should be swip_msstats
swip_msstats <- left_join(filt_results,msstats_prot,by="Protein")

# split into BioFraction's for ease of comprehension
# extract fraction from label
fraction <- regmatches(filt_results$Label,
		       regexpr("\\F[0-9]{1,2}",filt_results$Label))
filt_results <- tibble::add_column(filt_results,
				   "Fraction" = fraction,
				   .before="Protein")
results_list <- split(filt_results,fraction)

#length(results_list)
# [1] 7
#names(results_list)
# [1] "F10" "F4"  "F5"  "F6"  "F7"  "F8"  "F9"

# sort
results_list <- results_list[c("F4","F5","F6","F7","F8","F9","F10")]

# sort each df by pvalue and drop rows with issues or NA pvals
clean_results <- function(x) { # x = results_list[[1]]
	is_issue <- !is.na(x$issue)
	x <- x[!is_issue,] # drop rows with issue
	x$issue <- NULL
	x <- x[order(x$pvalue),]
	return(x)
}
final_results <- lapply(results_list, clean_results)

# save as excel worksheet
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(final_results, myfile)
message(paste("\nSaved",myfile))

# save as rda
myfile <- file.path(root,"rdata","swip_msstats.rda")
save(swip_msstats, file=myfile,version=2)
message(paste("\nSaved",myfile))

# any protein overlap?
sigprot_list <- lapply(results_list,function(x) {
	      x %>% filter(adj.pvalue<0.05) %>% 
		      dplyr::select(Protein) %>% 
		      unlist() %>% as.character() %>% unique() })

top_sigprots <- Reduce(intersect, sigprot_list)
idx <- match(top_sigprots, gene_map$uniprot)
message(c("\nProteins that are differentially abundant in all fractions: ",
	      paste(gene_map$symbol[idx],collapse=", ")))

# Examine the number of significant proteins
message(paste("\nNumber of differentially abundant (FDR < 0.05) proteins:"))
df <- sapply(results_list,function(x) sum(x$adj.pvalue<0.05,na.rm=TRUE))
knitr::kable(t(df))

## FIXME: flip sign of log2FC
## add percent WT
## Change Label to Contrast or Comparison
## add protein level data?
