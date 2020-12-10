#' lmerTestContrast
#'
#' @export lmerTestContrast
#'
#' @import lmerTest
#'
#' @importFrom dplyr %>%

lmerTestContrast <- function(fm, contrast,
                             df_prior = 0, s2_prior = 0) {

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
  # 	    setNames(nm = names(lme4::fixef(fm)))
  #   contrast[grepl(".*Mutant.*F4",names(contrast))] <- -1
  #   contrast[grepl(".*Control.*F4",names(contrast))] <- +1
  #
  #   # test the comparison!
  #   lmerTestContrast(fm, contrast)

  stopifnot(all(names(contrast) == names(lme4::fixef(fm))))

  pos_group <- names(contrast[contrast > 0])
  neg_group <- names(contrast[contrast < 0])
  comparison <- paste(pos_group, neg_group, sep = "-")

  # compute Satterthwaite degrees of freedom
  model_summary <- summary(fm, ddf = "Satterthwaite")

  # variance associated with mixed effects
  #mixef_var <- as.data.frame(lme4::VarCorr(fm, comp = "Variance"))

  # collect model's degrees of freedom and coefficients (beta)
  s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])[1]
  coeff <- model_summary$coeff[, "Estimate"] # == lme4::fixef(fm)

  # compute the unscaled covar matrix
  # all(unscaled_vcov * sigma^2 == vcov)
  unscaled_vcov <- fm@pp$unsc()

  # compute scaled variance-covariance matrix
  vcov <- as.matrix(model_summary$vcov) # == vcov(fm) == fm@vcov_beta

  # compute 
  se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance

  # extract asymtoptic var-covar matrix from fit model
  A <- fm@vcov_varpar

  # calculate gradient from gradient matrices
  # consider using lmerTest::calcSatterth(tv, L)
  g <- sapply(fm@Jac_list, function(gm) contrast %*% gm %*% contrast)

  # given gradient and asymptoptic var-covar, compute posterior df
  denom <- as.numeric(g %*% A %*% g)
  df_post <- 2 * se2^2 / denom + df_prior

  # calculate posterior s2
  s2 <- sigma(fm)^2
  s2_post <- (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)

  # compute variance given s2_post
  vcov_post <- unscaled_vcov * s2_post # sace as vcov
  variance <- as.numeric(contrast %*% vcov_post %*% contrast)

  # compute fold change and the t-statistic [lmerTest eq 11]
  FC <- (contrast %*% coeff)[, 1]
  t <- FC / sqrt(variance)

  # compute the p-value given t-statistic and posterior degrees of freedom
  p <- 2 * pt(-abs(t), df = df_post)

  # collect stats
  prot_stats <- data.frame(
    Contrast = comparison,
    log2FC = FC,
    percentControl = 2^FC,
    SE = sqrt(variance),
    Tstatistic = t,
    Pvalue = p,
    DF = df_post,
    S2 = s2,
    isSingular = lme4::isSingular(fm)
  )

  return(prot_stats)
} # EOF
