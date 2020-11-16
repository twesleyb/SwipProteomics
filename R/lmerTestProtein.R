#' lmerTestProtein
#' @import dplyr lmerTest
#' @export lmerTestProtein lmerTestContrast


lmerTestContrast <- function(fm, contrast, variance = FALSE,
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

  # stopifnot(sum(contrast[contrast<0],contrast[contrast>0])==0)

  pos_group <- names(contrast[contrast > 0])
  neg_group <- names(contrast[contrast < 0])
  comparison <- paste(pos_group, neg_group, sep = "-")

  # compute Satterthwaite degrees of freedom
  model_summary <- summary(fm, ddf = "Satterthwaite")

  # variance associated with mixed effects
  mixef_var <- as.data.frame(lme4::VarCorr(fm, comp = "Variance"))

  # collect model's degrees of freedom and coefficients (Beta)
  # s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])[1]
  coeff <- model_summary$coeff[, "Estimate"] # == lme4::fixef(fm)

  # NOTE: we just use the model summary's vcov object
  # compute the unscaled covar matrix with unsc() accessed within the model's
  # environment
  # unscaled_vcov <- fm@pp$unsc()
  # compute scaled variance-covariance matrix
  # all(unscaled_vcov * sigma^2 == vcov)

  # compute scaled variance-covariance matrix
  vcov <- as.matrix(model_summary$vcov) # == vcov(fm) == fm@vcov_beta

  # compute variance
  se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance

  # extract asymtoptic var-covar matrix from fit model
  A <- fm@vcov_varpar

  # calculate gradient from gradient matrices
  # consider using lmerTest::calcSatterth(tv, L)
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
  prot_stats <- data.frame(
    Contrast = comparison,
    log2FC = FC,
    percentControl = 2^FC,
    SE = sqrt(se2),
    Tstatistic = t,
    Pvalue = p,
    DF = df_post,
    isSingular = lme4::isSingular(fm)
  )

  if (variance) {
    # return the variance estimates of mixed effects
    mixef_var <- as.data.frame(lme4::VarCorr(fm, comp = "Variance"))
    attr(prot_stats, "variance") <- setNames(mixef_var$"vcov", nm = mixef_var$grp)
  }

  return(prot_stats)
} # EOF


lmerTestProtein <- function(protein, fx, msstats_prot, contrasts) {
  suppressPackageStartupMessages({
    library(lme4)
    library(dplyr)
    library(lmerTest)
  })

  getIndex <- function(namen, dm = lme4::fixef(fm)) {
    # a helper function to find column index
    idy <- which(grepl(namen, names(dm)))
  }

  # check input args
  stopifnot(inherits(fx, "formula"))
  stopifnot(inherits(protein, "character"))
  stopifnot(inherits(msstats_prot, "data.frame"))

  # input contrasts should be astats_df matrix, numeric vector, or list of such
  if (inherits(contrasts, "matrix")) {
    # rows of matrix are contrasts
    # rownames of matrix are contrast names (comparisons)
    contrast_list <- unlist(apply(contrasts, 1, list), recursive = FALSE)
    stopifnot(all(sapply(contrast_list, sum) == 0))
  } else if (inherits(contrasts, "numeric")) {
    # the sum of positive and neg coeff in contrast should be 0
    stopifnot(sum(contrasts[contrasts < 0], contrasts[contrasts > 0]) == 0)
    contrast_list <- list(contrasts)
  } else if (inherits(contrasts, "list")) {
    contrast_list <- contrasts
    # comparison names will be generated from pos and neg coefficients
  } else {
    stop("problem parsing input 'contrasts'.")
  }

  # subset the data
  idx <- msstats_prot$Protein == protein
  subdat <- msstats_prot[idx, ]
  if (any(is.na(subdat))) {
    # the data cannot contain missing values
    warning("The data cannot contain missing values.")
    return(NULL)
  }

  # fit the model
  fm <- lmerTest::lmer(fx, data = subdat)

  # evaluate statistical comparisons for all contrasts
  stats_list <- list()

  for (i in seq(contrast_list)) {

    # insure contrast matches names(fixef(fm))
    comparison <- contrast_list[[i]]
    contrast <- comparison[names(sort(sapply(names(comparison), getIndex)))]

    # assess contrast
    test_results <- lmerTestContrast(fm, contrast)

    # collect results
    test_results$Protein <- protein

    if (!is.null(attr(contrast_list, "names"))) {
      # use names of contrasts in contrast_list if they exist
      test_results$Contrast <- names(contrast_list)[i]
    }

    stats_list[[i]] <- test_results
  } # EOL for every comparison

  # compile results
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  # sort cols
  stats_df <- stats_df[, c(
    "Protein", "Contrast", "log2FC", "percentControl",
    "Tstatistic", "Pvalue", "SE", "DF", "isSingular"
  )]

  return(stats_df)
} # EOF
