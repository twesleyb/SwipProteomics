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

# load renv
renv::load(root)

# load projects functions in root/R
devtools::load_all()

# other imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
	library(MSstatsTMT)
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot

## munge sample to group mapping -----------------------------------------------

# munge groups for Control vs Mutant comparison
msstats_prot$BioReplicate <- gsub("\\.R[1,2,3]","",
				  as.character(msstats_prot$BioReplicate))
msstats_prot$BioReplicate <- as.factor(msstats_prot$BioReplicate)
msstats_prot$Condition <- gsub("\\.F[0-9]{1,2}$","",
			       as.character(msstats_prot$Condition))
msstats_prot$Condition <- as.factor(msstats_prot$Condition)

# annotate samples with biofraction?
#fraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
#condition <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
#msstats_prot$BioFraction <- factor(fraction,levels=c("F4","F5","F6","F7","F8","F9","F10"))
#msstats_prot$Condition <- factor(condition,levels="Mutant","Control")
#msstats_prot$BioReplicate <- msstats_prot$Run


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
load(file.path(root,"rdata","msstats_contrasts.rda"))
cm0 <- msstats_contrasts # intrafraction contrasts

# new contrast matrix for 'WT v Mutant' contrast adjusted for fraction
# differences
cm1  <- matrix(c(-1,1),nrow=1,ncol=2)
colnames(cm1) <- c("Control","Mutant")
rownames(cm1) <- "Mutant-Control"

proteins <- unique(as.character(msstats_prot$Protein))
prots <- sample(proteins,100)
subdat <- msstats_prot %>% filter(Protein %in% prots)

mstats_quant <- MSstats::groupComparison(cm1,subdat)

## models -------------------------------------------------

# Model fit by MSstatsTMT for intrafraction comparisons:
# 0. lmer: Abundance ~ (1|Run) + Condition
fx0 <- formula("Abundance ~ 1 + (1|Run) + Condition")

# Time-course (repeated measures) models fit by MSstats:
# 1. lmer: Abundance ~ Condition + (1|BioReplicate)
fx1 <- formula("Abundance ~ Condition + (1|BioReplicate)")

# 2. lmer: Abundance ~ Condition + (1|BioReplicate) + (1|Condition:BioReplicate)
fx2 <- formula("Abundance ~ Condition + (1|BioReplicate) + (1|Condition:BioReplicate)")


## fit models for swip --------------------------------------------------------

prot <- "Q3UMB9" # Swip

fit_list <- MSstatsTMT::fitLMER(fx0, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(contrast_list)

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

# get results list -- df for each comparison and compute padjust
results_list <- getResults(fit_list)

knitr::kable(bind_rows(results_list))
