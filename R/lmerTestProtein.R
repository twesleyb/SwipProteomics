#' lmerTestProtein
#' @import lmerTest
#' @export lmerTestContrast lmerTestProtein

lmerTestContrast <- function(fm, contrast, 
			     df_prior=0,s2_prior=0) {

    # define the comparison to be tested
    # FIXME: what if contrast is more complicated?
    pos_group <- names(contrast[contrast==1])
    neg_group <- names(contrast[contrast==-1])
    comparison <- paste(pos_group,neg_group,sep="-")

    # get model's summary and compute degrees of freedom
    model_summary <- summary(fm, ddf = "Satterthwaite")

    # collect model df and coefficients
    s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])
    coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)

    # extract sigma and theta
    sigma <- as.numeric(model_summary$sigma) # == stats::sigma(fm)
    theta <- as.numeric(model_summary$optinfo$val) # == lme4::getME(fm, "theta")

    # We need to compute the unscaled variance-covariance matrix!
    # Obscure, but here it is:
    #?vcov.merMod -- by inspection we find the calculation of scaled var-covar:
    #V <- sigm^2 * object@pp$unsc()
    # fm@pp$unsc is the unscaled variance covariance matrix!
    #vcov <- model_summary$vcov # == vcov(fm) == fm@vcov_beta
    unscaled_vcov <- fm@pp$unsc()

    #https://rdrr.io/cran/merDeriv/man/vcov.lmerMod.html
    # variance of the random effects
    #lme4::VarCorr(fm,sigma)

    #https://stackoverflow.com/questions/27113456/what-is-the-unscaled-variance-in-rs-linear-model-summary
    # unscaled: solve(t(X) %*% X)
    # scaled: solve(t(X) %*% X)*sigma^2 -- scaled by the variance sigma^2!

    # the variation not captured by fixed effects
    # attr(lme4::VarCorr(fm), "sc")^2

    # extract asymtoptic var-covar matrix from fit model
    # NOTE: done en route to se2 
    A <- fm@vcov_varpar

    # vcov_varpara numeric matrix holding the asymptotic variance-covariance 
    # matrix of the variance parameters (including sigma).
    # we compute the unmoderated statistics:
    # things that depend upon the contrast
    se2 <- as.numeric(contrast %*% unscaled_vcov %*% contrast) # == variance 

    # calculate posterior s2
    s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)

    # calculate gradient 
    g <- sapply(fm@Jac_list, function(gm) contrast %*% gm %*% contrast)

    # Jac_list - a list of gradient matrices (Jacobians) for the gradient of the
    # variance-covariance of beta with respect to the variance parameters, where
    # beta are the mean-value parameters available in fixef(object).
    # given gradient and asymptoptic var-covar
    
    # calculate denom
    denom <- as.numeric(t(g) %*% A %*% g)

    # compute posterior df
    # bayesian? posterior proporional to inverse of - σ² ?
    #df_post <- (2 *( se2/denom) + df_prior)
    df_post <- (se2 / denom) + df_prior

    # compute fold change and the t-statistic
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(se2) # 

    # compute the p-value given t-statistic and df.post
    p <- 2 * pt(-abs(t), df = df_post) 

    # collect stats
    prot_stats <- data.frame(Contrast=comparison,
			     log2FC=FC, 
			     percentControl=2^FC, 
			     Pvalue=p,
			     Tstatistic=t, 
			     SE=sqrt(se2), 
			     DF=df_post, 
			     isSingular=lme4::isSingular(fm))

    return(prot_stats)
} #EOF


lmerTestProtein <- function(protein, fx, msstats_prot, contrasts) {

  # input contrasts should be a numeric vector, or list of such 
  if (inherits(contrasts,"numeric")) {
  	contrasts <- list(contrasts)
  } else {
	stopifnot(inherits(contrasts,"list"))
  }

  # subset the data
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
	# the data cannot contain missing values
  	warning("The data cannot contain missing values.")
        return(NULL)
  }

  # fit the model
  fm <- lmerTest::lmer(fx, data=subdat)

  # evaluate statistical comparisons for all contrasts
  stats_list <- list()
  for (i in seq(contrasts)) {
	  contrast <- contrasts[[i]]
	  df <- lmerTestContrast(fm, contrast)
	  df$Protein <- protein # add protein annotation!
	  stats_list[[i]] <- df
  }

  # compile results
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL
  # sort cols
  stats_df[,c("Protein","Contrast","log2FC","percentControl",
	      "Pvalue","Tstatistic","SE", "DF", "isSingular")]

  return(stats_df)
} #EOF
