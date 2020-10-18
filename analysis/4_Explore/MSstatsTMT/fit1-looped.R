#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

# input:
root <- "~/projects/SwipProteomics"
# * msstats_prot

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root,quiet=TRUE) }

# imports
suppressPackageStartupMessages({
  library(dplyr)
})


## Functions ------------------------------------------------------------------

fit1 <- function(protein) {
  #############################################################################
  ## Q1. Which model is the correct model?
  ## We are interested in differences between levels of Genotype
  ## adjusted for 'BioFraction'.
  ## [1] all covariates
  # boundary fit -- problem with combination of mixture and subject
  #fx <- formula("Abundance ~ BioFraction + (1|Mixture) + (1|Subject) + Genotype") 
  ## [2] Only mixed effect of Subject:
  fx <- formula("Abundance ~ BioFraction + (1|Subject) + Genotype")
  #r.squaredGLMM.merMod(fm)
  #    R2m       R2c
  # [1,] 0.9162647 0.9663459 <-- better fit(?) (more var explained)
  # R2c is total variance explained
  # R2m is the variance explained by fixed effects
  # remaineder = ~0.05 is the variance explained by mixed effects.
  ## [3] Only mixed effect of Mixture:
  #fx <- formula("Abundance ~ BioFraction + (1|Mixture) + Genotype")
  #r.squaredGLMM.merMod(fm)
  #    R2m       R2c
  #[1,] 0.9222618 0.9355818  <-- R2c = total variance explained 
  #############################################################################
  # the data cannot contain missing values
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
  	warning("The data cannot contain missing values.")
        return(NULL)
  }
  # fit the model
  opts <- lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "stop", tol=1e-4))
  fm <- try(lmerTest::lmer(fx, data=subdat, control = opts),silent=TRUE)
  if (inherits(fm,"try-error")) { 
	  warning(protein ," results in a singular fit, returning NULL.")
	  return(NULL) 
  }
  #############################################################################
  ## Q3. please help me better understand output of:
  #summary(fm)
  #############################################################################
  #############################################################################
  ## Q4. Goodness of fit. How to asses?
  #qqnorm(resid(fm))
  #qqline(resid(fm))
  #############################################################################
  # compute some model statistics, store in list rho
  rho <- list()
  rho$protein <- protein
  rho$formula <- fx
  rho$model <- fm
  rho$data <- subdat
  # here we compute the unmoderated statistics
  rho$df.prior <- 0
  rho$s2.prior <- 0
  # check if singular (boundary) fit
  rho$isSingular <- lme4::isSingular(fm)
  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$thopt <- lme4::getME(fm, "theta")
  # calculate degrees of freedom and sigma^2:
  av <- anova(fm)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"
  # calcuate symtoptic var-covar matrix
  rho$A <- MSstatsTMT::calcApvarcovar(fx, subdat, rho$thopt,rho$sigma) 
  # we store the proteins statistics in a list
  stats_list <- list()
  # calculate posterior s2
  s2.post <- MSstatsTMT::calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)
  # compute variance-covariance matrix
  vss <- MSstatsTMT::vcovLThetaLM(fx,subdat)
  #############################################################################
  ## Q5. How to define contrast matrix?
  vec = rho$coeff
  vec[] <- 0
  #vec[1] <- -1
  vec[length(vec)] <- 1
  contrast_matrix <- vec
  #############################################################################
  varcor <- vss(t(contrast_matrix), thpars=c(rho$thopt, rho$sigma))
  vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
  se2 <- as.numeric(t(contrast_matrix) %*% as.matrix(vcov) %*% contrast_matrix)
  # calculate variance
  vcov.post <- varcor$unscaled.varcor * s2.post
  variance <- as.numeric(t(contrast_matrix) %*% as.matrix(vcov.post) %*% contrast_matrix)
  # calculate degrees of freedom
  # given params theta and sigma from lmer and the Acovar
  # NOTE: careful formatting can break things
  g <- MSstatsTMT::mygradient(function(x) vss(t(contrast_matrix), x)$varcor, c(rho$thopt, rho$sigma))
  denom <- t(g) %*% rho$A %*% g
  # compute df.posterior
  df.post <- 2 * (se2)^2 / denom + rho$df.prior # df.post
  ## Q5. FC seems inflated?
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% rho$coeff)[, 1] # coeff
  t <- FC / sqrt(variance) # the Fold change and the variance
  # compute the p-value
  p <- 2 * pt(-abs(t), df = df.post) # t-statistic and df.post
  # compile results
  rho$stats <- data.frame(protein=protein,contrast="Mutant-Control",
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
  		 Tstatistic=t, SE=sqrt(variance), DF=df.post, isSingular=rho$isSingular)
  return(rho)
} #EOF

## load the data --------------------------------------------------------------

# load msstats preprocessed protein data from SwipProteomics in root/data
#devtools::load_all()
load(file.path(root,"data","swip.rda"))
load(file.path(root,"data","gene_map.rda"))
load(file.path(root,"data","msstats_prot.rda"))

washc = gene_map$uniprot[grep("Washc*",gene_map$symbol)]

# Munge sample annotations:
# * create Genotype column
# * create BioFraction column
genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
subject <- interaction(msstats_prot$Mixture,genotype)
msstats_prot$Genotype <- as.factor(genotype)
msstats_prot$BioFraction <- biofraction
msstats_prot$Subject <- subject

results <- list()
#proteins = c(washc, sample(unique(as.character(msstats_prot$Protein)),100))
proteins = unique(as.character(msstats_prot$Protein))
# loop to fit protein-wise models and perform statistical comparisons
pbar <- txtProgressBar(max=length(proteins),style=3)
for (protein in proteins) {

	results[[protein]] <- fit1(protein)

	setTxtProgressBar(pbar, value = match(protein,proteins))
}
close(pbar)

# error:
#"P00375"
#Error in t(g) %*% rho$A : 
#requires numeric/complex matrix/vector arguments


# collect the results and perform Padjust
results_df <- bind_rows(sapply(results,"[","stats"))
results_df$FDR <- p.adjust(results_df$Pvalue,"BH")

# check
results_df %>% arrange(Pvalue) %>% knitr::kable()

sum(results_df$FDR< 0.1)
