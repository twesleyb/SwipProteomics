 #!/usr/bin/env Rscript 

# description: Working through MSstatsTMT protein-level models and 
# statistical comparisons

# input:
# * preprocessed protein-level data from PDtoMSstatsTMTFormat()
# comp - all intrafraction comparisons in the form "Mutant.F4-Control.F4"
# myfile <- file.path(root,"rdata","msstats_prot.rda")


## prepare environment --------------------------------------------------------

# project root dir and dir containing MSstatsTMT/R source code
root ="~/projects/SwipProteomics"
funcdir = "~/projects/SwipProteomics/src/MSstatsTMT"
#^ these are core MSstats internal functions used by MSstatsTMT_wrapper functions

# load renv
renv::load(root)

# load projects functions in root/R
# includes MSstatsTMT_wrappers.R
devtools::load_all()

# load MSstatsTMT's guts
load_fun(funcdir)

# other imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot


## begin protein-level modeling -------------------------------------------------

# for intra-fraction comparisons MSstatsTMT fits the model:
# FIXME: change run to batch or experiment... BATCH
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")
message("Fitting protein-wise mixed-effects linear model of the form:\n\t",fx)

# fit model for Swip
fit_list <- fitLMER(fx,msstats_prot,protein="Q3UMB9")

length(fit_list)
names(fit_list)
str(fit_list[[1]])

## build a contrast_matrix ----------------------------------------------------

# define all intrafraction comparisons
comp <- paste(paste("Mutant",paste0("F",seq(4,10)),sep="."),
	      paste("Control",paste0("F",seq(4,10)), sep="."), sep="-")

head(comp,3)

# create a contrast matrix
contrast_matrix <- getContrasts(comp, groups=levels(msstats_prot$Condition))

dim(contrast_matrix)
rownames(contrast_matrix)[1]
knitr::kable(contrast_matrix[1,])


## test contrasts -------------------------------------------------------------

# for each protein compare conditions declared in contrast_matrix
fit_list <- testContrasts(fit_list, contrast_matrix)

names(fit_list)

fit_list[[1]]
