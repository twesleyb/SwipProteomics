#!/usr/bin/env Rscript

# load renv
root = "~/projects/SwipProteomics"
renv::load(root)

# import
suppressPackageStartupMessages({
  library(dplyr)
  library(MSstatsTMT)
})

devtools::load_all() # and internal funcs
# FIXME: what internals are required?

# load msstats preprocessed protein data saved in root/rdata
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile) # == msstats_prot


## functions ------------------------------------------------------------------

getRho <- function(fm) {
  #' getRho 
  #' utilized by fitLMER to calculate some statistics from the fit model
  rho <- list()
  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$thopt <- lme4::getME(fm, "theta")
  # calculate degrees of freedom and sigma^2:
  av <- anova(fm)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"
  # calcuate symtoptic var-covar matrix
  # NOTE: utlilizes MSstatsTMT internal function to do calculation
  rho$model <- fm
  # FIXME need to put calcApvar here!
  rho$A <- .calcApvar(rho) #
  return(rho)
} # EOF


calcPosterior <- function(s2, s2_df, s2.prior = 0, df.prior = 0) {
  # A function to compute s2 posterior
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}


## START ----------------------------------------------------------------------

protein = "Q3UMB9"
moderated = FALSE

# If protein == NULL then all proteins will be fit
if (!is.null(protein)) {
  check <- is.character(protein) &
  protein %in% msstats_prot$Protein
  stopifnot(check)
  # the specified proteins will be fit
  proteins <- protein
  } else if (is.null(protein)) {
    # all proteins in msstats_prot will be fit
    proteins <- unique(as.character(msstats_prot$Protein))
  } else {
    stop("problem parsing input 'protein'")
} #EIS

## the data cannot contain missing values
if (any(is.na(msstats_prot %>% filter(Protein %in% proteins)))) {
# FIXME: shouldn't missing values have been imputed by MSstats prev?
to_drop <- unique(msstats_prot$Protein[is.na(msstats_prot$Abundance)])
msstats_prot <- msstats_prot[!is.na(msstats_prot$Abundance), ]
warning("Missing values are not tolerated in input 'data'.",
    "\n",length(to_drop)," proteins were removed.")
}

# subset the data
protein <- sample(unique(as.character(msstats_prot$Protein)),1)
subdat <- msstats_prot %>% filter(Protein == protein)
# Munge sample annotations
subdat$Genotype <- as.factor(sapply(strsplit(as.character(subdat$Condition),"\\."),"[",1))
subdat$BioFraction <- as.factor(sapply(strsplit(as.character(subdat$Condition),"\\."),"[",2))

# msstats fits the model:
# ABUNDANCE ~ GROUP + (1|SUBJECT) | Group is ~Condition and SUBJECT ~ Fraction

# fit the model
# SOMETHING SEEMS WRONG ABOUT THIS MODEL. 
fx <- formula("Abundance ~ + (1|BioFraction) + Genotype")
message("lmer: ",fx)
fm <- lmerTest::lmer(fx, data=subdat)

fm

rho <- getRho(fm)
rho$data <- subdat
rho[["isSingular"]] <- lme4::isSingular(fm)
# contrast matrix for comparison to control
cm <- rho$coeff
cm[1] <- -1
cm[2] <- 1
# compute prior df and s2
if (!moderated) {
  # the unmoderated t-statistic:
  df.prior <- 0
  s2.prior <- 0
  } else if (moderated) {
    idx <- sapply(fit_list, "[", "s2_df") != 0
    eb_input_s2 <- as.numeric(sapply(fit_list, "[", "s2")[idx])
    eb_input_df <- as.numeric(sapply(fit_list, "[", "s2_df")[idx])
    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    if (is.infinite(eb_fit$df.prior)) {
      df.prior <- 0
      s2.prior <- 0
    } else {
      df.prior <- eb_fit$df.prior
      s2.prior <- eb_fit$var.prior
  }
} # EIS fin calc of prior df and s2
# we store the proteins statistics in a list
stats_list <- list()
stats_list$Protein <- protein
# compute s2 posterior
s2.post <- calcPosterior(s2 = rho$s2, s2_df = rho$s2_df)
# compuate variance-covariance matrix
vss <- .vcovLThetaL(fm)
varcor <- vss(t(cm), c(rho$thopt, rho$sigma))
vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
# calculate variance
vcov.post <- varcor$unscaled.varcor * s2.post
variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
# calculate degrees of freedom
g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
denom <- t(g) %*% rho$A %*% g
df.post <- 2 * (se2)^2 / denom + df.prior
# compute fold change and the t-statistic
FC <- (cm %*% rho$coeff)[, 1]
t <- FC / sqrt(variance)
# compute the p-value
p <- 2 * pt(-abs(t), df = df.post)
df <- data.frame("Protein"=protein,
		 "P-value"=as.numeric(p),
		 "Log2 FC" = as.numeric(FC),
		 "t-statistic" = as.numeric(t),
		 "variance"=as.numeric(variance))
knitr::kable(df)

print(summary(fm))
