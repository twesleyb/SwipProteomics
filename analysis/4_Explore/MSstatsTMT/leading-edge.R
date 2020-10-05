 #!/usr/bin/env Rscript 

# description: Working through MSstatsTMT protein-level models and 
# statistical comparisons

# input:
# * preprocessed protein-level data from PDtoMSstatsTMTFormat()
# comp - all intrafraction comparisons in the form "Mutant.F4-Control.F4"
# myfile <- file.path(root,"rdata","msstats_prot.rda")


## prepare environment --------------------------------------------------------

# project root dir:
root ="~/projects/SwipProteomics"

# MSstats source directory:
funcdir = "~/projects/SwipProteomics/src/MSstatsTMT"
#^ these are core MSstats internal functions used by MSstatsTMT_wrapper functions

# load renv
renv::load(root)

# load projects functions in root/R
# includes MSstatsTMT_wrappers.R utilized herein
devtools::load_all()

# load MSstatsTMT's guts
load_fun(funcdir)

# other imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
	#library(doParallel) # for parallel processing
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot

# munge for Control vs Mutant comparison
msstats_prot$BioReplicate <- gsub("\\.R[1,2,3]","",
				  as.character(msstats_prot$BioReplicate))
msstats_prot$BioReplicate <- as.factor(msstats_prot$BioReplicate)
msstats_prot$Condition <- gsub("\\.F[0-9]{1,2}$","",
			       as.character(msstats_prot$Condition))
msstats_prot$Condition <- as.factor(msstats_prot$Condition)


# annotate samples with biofraction
#fraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
#condition <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
#msstats_prot$BioFraction <- fraction
#msstats_prot$Condition <- condition


## build a contrast_matrix ----------------------------------------------------

# pepare a contrast matrix
contrast_matrix <- matrix(c(-1,1),nrow=1,ncol=2)
colnames(contrast_matrix) <- c("Control","Mutant")
rownames(contrast_matrix) <- "Mutant-Control"


## begin protein-level modeling -------------------------------------------------

# for time-course aka repeated measures designs MSstatsTMT fits the model:
# FIXME: change run to batch or experiment... BATCH
fx <- formula("Abundance ~ 1 + (1|BioReplicate) + Condition")
message("Fitting protein-wise mixed-effects linear model of the form:\n\t",fx)

# FIXME: timit! parallize!
prots <- c("Q3UMB9", sample(unique(as.character(msstats_prot$Protein)),9))
fit_list <- MSstatsTMT::fitLMER(fx,msstats_prot,protein=prots)

# FIXME: how to catch warnings from 
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  : 
#Model failed to converge with max|grad| = 0.00223701 (tol = 0.002, component 1)
#fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, protein=prot)


## test contrasts -------------------------------------------------------------

# for each protein compare conditions declared in contrast_matrix
# FIXME: if run twice then errors
# Error in contrast.single[contrast.single > 0] <- temp :
# NAs are not allowed in subscripted assignments 
# FIXME: contrast.single problem when changing design
fit_list <- MSstatsTMT::testContrasts(fit_list, contrast_matrix, moderated = TRUE)

# get results list -- df for each comparison and compute padjust
results_list <- adjustPvalues(fit_list)

knitr::kable(bind_rows(results_list))
