#!/usr/bin/env Rscript

# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(swip)
data(msstats_prot)

# other imports
library(dplyr)
library(data.table)
library(doParallel)

## functions 

lmerFit <- function(fx,msstats_prot,protein) {
	# fit the model and catches errors
	fm <- suppressMessages({
		try(lmerTest::lmer(fx,msstats_prot %>% 
				 filter(Protein == protein)),silent=TRUE)
	})
	if (inherits(fm,"try-error")) {
		return(NULL)
	} else {
		return(fm)
	}
} # EOF


## fit model to every protein -------------------------------------------------

# the formula to be fit:
fx <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)")

# loop to fit model to every protein, extract coeff (fixed-effects)
proteins <- unique(as.character(msstats_prot$Protein))

n <- parallel::detectCores() -1 
doParallel::registerDoParallel(n)

# fit the models in parallel
fit_list <- foreach(protein = proteins) %dopar% {
	lmerFit(fx,msstats_prot,protein)
}
names(fit_list) <- proteins

message("Fit ", length(fit_list), " models.")


## clean-up fit list ----------------------------------------------------------
# remove null
idx <- !sapply(fit_list,is.null)
filt_list <- fit_list[idx]
message("Removed ", sum(!idx), " NULL models.")

# get fixed effects (coefficients)
fixef_list <- lapply(filt_list,lme4::fixef)

# drop incomplete models
idx <- sapply(fixef_list,function(x) length(x) == 14 & is.numeric(x))
fixef_list <- fixef_list[idx]
fm_list <- filt_list[idx]
message("Removed ", sum(!idx), " incomplete (rank-deficient) models.")


## gof ------------------------------------------------------------------------

# calculate goodness of fit
df <- as.data.table(t(sapply(fm_list, r.squaredGLMM.merMod)),
		    keep.rownames="Protein")
colnames(df)[colnames(df)=="V1"] <- "R2m" # fixed effects
colnames(df)[colnames(df)=="V2"] <- "R2c" # total

# remove proteins with poor fits 
# (percent var explained by fixef Geno:BioF < 0.5)
out <- df %>% filter(R2m<0.5) %>% select(Protein) %>% unlist() %>% unique()
keep <- names(fm_list)[names(fm_list) %notin% out]
message("Removed ", length(out), " models with poor fit.")


## create covariation networks ------------------------------------------------

# create data matrix
dm <- do.call(rbind,fixef_list[keep])

# calculate pearson correlation
adjm <- cor(t(dm))

# perform network enhancement
ne_adjm <- neten::neten(adjm)


## save the data --------------------------------------------------------------

# save adjm
myfile <- file.path(root,"rdata","adjm.rda")
save(adjm,file=myfile,version=2)

# save  ne adjm as rda
myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm,file=myfile,version=2)

# save ne adjm as csv for Leidenalg
myfile <- file.path(root,"rdata", "ne_adjm.csv")
ne_adjm %>% as.data.table(keep.rownames="Protein") %>% fwrite(myfile)
