#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(swip)
data(msstats_prot)

library(dplyr)
library(data.table)
library(doParallel)

lmerFit <- function(fx,msstats_prot,protein) {
	# just fits the model and catches failures
	fm <- try(lmerTest::lmer(fx,msstats_prot %>% 
				 filter(Protein == protein)),silent=TRUE)
	if (inherits(fm,"try-error")) {
		return(NULL)
	} else {
		return(fm)
	}
} # EOF



# loop to fit all models, extract coeff (fixed-effects)

proteins <- unique(as.character(msstats_prot$Protein))

fx <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)")

n <- parallel::detectCores() -1 
doParallel::registerDoParallel(n)

# fit the models in parallel
fit_list <- foreach(protein = proteins) %dopar% {
	lmerFit(fx,msstats_prot,protein)
}
names(fit_list) <- proteins


# remove null
idx <- which(!sapply(fit_list,is.null))
filt_list <- fit_list[idx]

# get fixed effects (coefficients)
fixef_list <- lapply(filt_list,lme4::fixef)

# drop incomplete models
idx <- which(sapply(fixef_list,function(x) length(x) == 14 & is.numeric(x)))
fixef_list <- fixef_list[idx]
fm_list <- filt_list[idx]

# calculate goodness of fit
df <- as.data.table(t(sapply(fm_list, r.squaredGLMM.merMod)),
		    keep.rownames="Protein")
colnames(df)[colnames(df)=="V1"] <- "R2m" # fixed effects
colnames(df)[colnames(df)=="V2"] <- "R2c" # total

out <- df %>% filter(R2m<0.5) %>% select(Protein) %>% unlist() %>% unique()
keep <- names(fm_list)[names(fm_list) %notin% out]

# create data matrix
dm <- do.call(rbind,fixef_list[keep])

# calculate pearson correlation
adjm <- cor(t(dm))

# perform network enhancement
ne_adjm <- neten::neten(adjm)

# save
myfile <- file.path(root,"rdata", "ne_adjm.csv")
ne_adjm %>% as.data.table(keep.rownames="Protein") %>% fwrite(myfile)

myfile <- file.path(root,"rdata","adjm.rda")
save(adjm,file=myfile,version=2)

save(ne_adjm,file=myfile,version=2)
myfile <- file.path(root,"rdata","ne_adjm.rda")
