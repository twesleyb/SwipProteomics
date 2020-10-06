#!/usr/bin/env Rscript 

# title:
# author: twab
# description: Working through MSstatsTMT protein-level models and 
#   statistical comparisons


## prepare environment --------------------------------------------------------

# project root dir:
root <- "~/projects/SwipProteomics"

# load renv
renv::load(root)

# other imports
suppressPackageStartupMessages({
	library(dplyr) 
	library(MSstatsTMT) # my fork
})

# local imports
devtools::load_all() 


## load the data ---------------------------------------------------------------

# load saved contrast matrix
data(msstats_contrasts)

# load msstats preprocessed protein data saved in root/rdata
data(msstats_prot)

# munge - clarify covariate names
mixture <- gsub("_1","",as.character(msstats_prot$Run))
genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)

msstats_prot$Mixture <- as.factor(mixture)
msstats_prot$Genotype <- as.factor(genotype)
msstats_prot$BioFraction <- as.factor(biofraction)


## fit model for swip --------------------------------------------------------

swip <- "Q3UMB9"

# fit the model for intrafraction comparisons
fx0 <- Abundance ~ 1 + (1|Mixture) + Condition
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

message("lmer:",fx0)
print(summary(fm0))

# fit the model for WT-Mutant comparisons
fx1 <- Abundance ~ 0 + (1|BioFraction) + Genotype
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein == swip))

message("lmer:",fx1)
print(summary(fm1))

# The statistical analysis for Swip, fx0
# NOTE: mstats_contrasts is all pairwise comparisons
fx0 <- Abundance ~ (1|Run) + Condition # This is how MSstatsTMT expects it to be
fit_list <- MSstatsTMT::fitLMER(fx0, msstats_prot, swip)
#fit_list <- MSstatsTMT::testContrasts(fit_list, msstats_contrasts)
#results <- getResults(fit_list)
#df <- bind_rows(results)
#knitr::kable(df)

# The statistical analysis for Swip, fx1
cm1 <- lme4::fixef(fm1)
cm1[1] <- -1 # Control
cm1[2] <- +1 # Mutant
fit_list <- MSstatsTMT::fitLMER(fx1, msstats_prot, swip)
fit_list <- MSstatsTMT::testContrasts(fit_list, cm1)
results <- getResults(fit_list)

#df <- bind_rows(results)
#knitr::kable(df)


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
data(msstats_contrasts)

cm <- lme4::coeff(fm0)
cm[1] <- -1
cm[2] <- 1
rho$contrasts <- cm


## fit model for swip --------------------------------------------------------

#fx0 <- formula("Abundance ~ 1 + (1|Mixture) + Condition")
fx1 <- formula("Abundance ~ 0 + (1|BioFraction) + Genotype")

fit <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(fit_list)

df <- bind_rows(results)

knitr::kable(df)

print(summary(fit))

names(fit_list[[1]])

fit_list[[1]]$coeff

fit_list[[1]]$sigma

fit_list[[1]]$s2_df

fit_list[[1]]$s2


fm <- lmerTest::lmer(fx, data=subdat)

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
# NOTE: utlilizes MSstatsTMT internal function .calcApvar to do calculation
rho$A <- MSstatsTMT::calcApvar(rho) 

knitr::kable(as.data.frame(rho$A))

# define contrast matrix for comparison to control
cm <- rho$coeff
cm[1] <- -1
cm[2] <- 1
rho$contrasts <- cm

message("Contrast:")
knitr::kable(t(cm))

# compute prior df and s2
#if (!moderated) {
#  # the unmoderated t-statistic:
#  df.prior <- 0
#  s2.prior <- 0
#  } else if (moderated) {
#    idx <- sapply(fit_list, "[", "s2_df") != 0
#    eb_input_s2 <- as.numeric(sapply(fit_list, "[", "s2")[idx])
#    eb_input_df <- as.numeric(sapply(fit_list, "[", "s2_df")[idx])
#    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
#    if (is.infinite(eb_fit$df.prior)) {
#      df.prior <- 0
#      s2.prior <- 0
#    } else {
#      df.prior <- eb_fit$df.prior
#      s2.prior <- eb_fit$var.prior
#  }
#} # EIS fin calc of prior df and s2

# we store the proteins statistics in a list
stats_list <- list()

# compute s2 posterior
calcPosterior <- function(s2, s2_df, s2.prior, df.prior) {
  # A function to compute s2 posterior
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}

s2.post <- calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)

# compuate variance-covariance matrix
vss <- vcovLThetaL(fm)
varcor <- vss(t(cm), c(rho$thopt, rho$sigma))
vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
se2 <- as.numeric(t(cm) %*% as.matrix(vcov) %*% cm)

# calculate variance
vcov.post <- varcor$unscaled.varcor * s2.post
variance <- as.numeric(t(cm) %*% as.matrix(vcov.post) %*% cm)

# calculate degrees of freedom
# given params theta and sigma from lmer and the Acovar
g <- mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
denom <- t(g) %*% rho$A %*% g

# compute df.posterior
df.post <- 2 * (se2)^2 / denom + rho$df.prior

# compute fold change and the t-statistic
FC <- (cm %*% rho$coeff)[, 1]
t <- FC / sqrt(variance)

# compute the p-value
p <- 2 * pt(-abs(t), df = df.post)

df <- data.frame("Protein"=protein,
		 "P-value"=as.numeric(p),
		 "Log2 FC" = as.numeric(FC),
		 "t-statistic" = as.numeric(t),
#!/usr/bin/env Rscript 

# title:
# author: twab
# description: Working through MSstatsTMT protein-level models and 
#   statistical comparisons


## prepare environment --------------------------------------------------------

# project root dir:
root <- "~/projects/SwipProteomics"

# load renv
renv::load(root)

# other imports
suppressPackageStartupMessages({
	library(dplyr) 
	library(MSstatsTMT) # my fork
})


## load the data ---------------------------------------------------------------

# load msstats preprocessed protein data saved in root/rdata
devtools::load_all() 
data(msstats_prot)

# munge - clarify covariate names
mixture <- gsub("_1","",as.character(msstats_prot$Run))
msstats_prot$Mixture <- mixture


## build a contrast_matrix ----------------------------------------------------

# load saved contrast matrix
data(msstats_contrasts)
cm0 <- msstats_contrasts 

# all intrafraction contrasts in the format MSstatsTMT expects
# example, a single row/contrast:
knitr::kable(cm0[1,])


## fit model for swip --------------------------------------------------------

prot <- "Q3UMB9"

fx <- formula("Abundance ~ 1 + (1|Mixture) + Condition")

fit <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == prot))

fit_list <- MSstatsTMT::fitLMER(fx, msstats_prot, prot)

fit_list <- MSstatsTMT::testContrasts(fit_list, cm0)

results <- getResults(fit_list)

df <- bind_rows(results)

knitr::kable(df)

print(summary(fit))

names(fit_list[[1]])

fit_list[[1]]$coeff

fit_list[[1]]$sigma

fit_list[[1]]$s2_df

fit_list[[1]]$s2

fit_list[[1]]$A
