#!/usr/bin/env Rscript

# Two steps:
# (1) given tidy, normalized protein data, fit the appropriate LMM (lmer)
# (2) given the model, assess a contrast (lmerTestContrast)

## functions ------------------------------------------------------------------

#' lmerTestContrast
#' @export lmerTestContrast 

lmerTestContrast <- function(fm, contrast,
			     df_prior=0,s2_prior=0) {
#    # example:
#    library(dplyr)
#    library(SwipProteomics)
#
#    data(msstats_prot)
#    data(swip) # "Q3UMB9"
#
#    fx <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)")
#    fm <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == swip))
#
#    # create a contrast!
#    geno <- msstats_prot$Genotype
#    biof <- msstats_prot$BioFraction
#    conditions <- levels(interaction(geno, biof))
#    contrast <- vector("numeric", length(conditions)) %>%
#	    setNames(nm = names(lme4::fixef(fm)))
#   contrast[grepl(".*Mutant.*F4",names(contrast))] <- -1
#   contrast[grepl(".*Control.*F4",names(contrast))] <- +1
#
#   # test the comparison!
#   lmerTestContrast(fm, contrast) 

    # comparison should be a numeric and sum of vector == 0
    stopifnot(inherits(contrast,"numeric"))
    stopifnot(sum(contrast[contrast<0],contrast[contrast>0])==0)
    pos_group <- names(contrast[contrast>0])
    neg_group <- names(contrast[contrast<0])
    comparison <- paste(pos_group,neg_group,sep="-")

    # compute Satterthwaite degrees of freedom
    model_summary <- summary(fm, ddf = "Satterthwaite")

    # variance associated with mixed effects
    mixef_var <- as.data.frame(lme4::VarCorr(fm,comp="Variance"))

    # collect model's coefficients (beta)
    coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)

    # compute scaled variance-covariance matrix
    vcov <- as.matrix(model_summary$vcov) # == vcov(fm) == fm@vcov_beta

    # compute variance
    se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance 

    # extract asymtoptic var-covar matrix from fit model
    A <- fm@vcov_varpar

    # calculate gradient from gradient matrices
    g <- sapply(fm@Jac_list, function(gm) contrast %*% gm %*% contrast)

    # given gradient and asymptoptic var-covar, compute posterior df
    denom <- as.numeric(g %*% A %*% g)
    df_post <- 2 * (se2^2 / denom) + df_prior 

    # compute fold change and the t-statistic [lmerTest eq 11]
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(se2) 

    # compute the p-value given t-statistic and posterior degrees of freedom
    p <- 2 * pt(-abs(t), df = df_post) 

    # collect stats
    prot_stats <- data.frame(Contrast=comparison,
			     log2FC=FC, 
			     percentControl=2^FC, 
			     SE=sqrt(se2), 
			     Tstatistic=t, 
			     Pvalue=p,
			     DF=df_post, 
			     isSingular=lme4::isSingular(fm))
    return(prot_stats)
} #EOF


#' getContrast
#' @export getContrast 

getContrast <- function(fm, pos_coef, neg_coef){

  # create a contrast!
  stopifnot(inherits(fm,"lmerModLmerTest"))
  stopifnot(is.character(pos_coef))
  stopifnot(is.character(neg_coef))

  contrast <- lme4::fixef(fm)

  contrast[] <- 0
  neg_index <- grepl(neg_coef,names(contrast))
  pos_index <- grepl(pos_coef,names(contrast))
  contrast[neg_index] <- -1/sum(neg_index)
  contrast[pos_index] <- +1/sum(pos_index)

  stopifnot(any(contrast > 0))
  stopifnot(any(contrast < 0))
  stopifnot(sum(contrast)==0)

  return(contrast)
}

## main -----------------------------------------------------------------------

# for mixed models and stats
library(lme4)
library(lmerTest)

# for manipulating the data
library(dplyr)
library(data.table)


# data 
data(swip) 
data(washc_prots)
data(msstats_prot)


# a mixed model formula or function (fx) for protein-level comparisons
fx0 <- "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)"

# fit the model (fm) to a subset of the data for a single protein
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% subset(Protein == swip))

# see more info about the model, compute Satterthwaite degrees of freedom
summary(fm0, ddf="Satterthwaite")


## examine two types of contrast:

# intra-BioFraction:
L1 <- getContrast(fm0,"GenotypeMutant:BioFractionF4","GenotypeControl:BioFractionF4")

# Mutant-Control:
L2 <- getContrast(fm0,"GenotypeMutant","GenotypeControl")

res1 <- lmerTestContrast(fm0,L1)

res2 <- lmerTestContrast(fm0,L2)

rbind(res1,res2) %>% knitr::kable()
