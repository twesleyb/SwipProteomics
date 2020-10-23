#' lmerTestProtein
#' @import dplyr lmerTest
#' @export lmerTestContrast lmerTestProtein


lmerTestContrast <- function(fm, contrast, 
			     df_prior=0,s2_prior=0) {

    # define the comparison to be tested
    stopifnot(sum(as.numeric(contrast))==0)
    pos_group <- names(contrast[contrast==1])
    neg_group <- names(contrast[contrast==-1])
    comparison <- paste(pos_group,neg_group,sep="-")

    # compute Satterthwaite degrees of freedom
    model_summary <- summary(fm, ddf = "Satterthwaite")

    # collect model's df and coefficients
    s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])
    coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)

    # extract sigma and theta
    sigma <- as.numeric(model_summary$sigma) # == stats::sigma(fm)
    theta <- as.numeric(model_summary$optinfo$val) # == lme4::getME(fm, "theta")

    # calculate posterior s2
    s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)

    # compute unscaled variance-covariance matrix
    #?vcov.merMod -- by inspection we find the calculation of scaled var-covar:
    # >>> V <- sigm^2 * object@pp$unsc()
    # We can compute the unscaled covar matrix with unsc() accessed within the 
    # models environment:
    unscaled_vcov <- fm@pp$unsc()

    # compute scaled variance-covariance matrix
    vcov <- as.matrix(model_summary$vcov) # == vcov(fm) == fm@vcov_beta

    # compute variance
    vcov_post <- unscaled_vcov * s2_post
    variance <- as.numeric(contrast %*% vcov_post %*% contrast)

    # extract asymtoptic var-covar matrix from fit model
    A <- fm@vcov_varpar

    # calculate gradient 
    g <- sapply(fm@Jac_list, function(gm) contrast %*% gm %*% contrast)

    # Jac_list - a list of gradient matrices (Jacobians) for the gradient of the
    # variance-covariance of beta with respect to the variance parameters, where
    # beta are the mean-value parameters available in fixef(object).
    # given gradient and asymptoptic var-covar
    
    # given gradient and asy var-covar, compute posterior df
    denom <- as.numeric(t(g) %*% A %*% g)
    se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance 
    df_post <- (se2 / denom) + df_prior

    # compute fold change and the t-statistic
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(variance) # 

    # compute the p-value given t-statistic and df.post
    p <- 2 * pt(-abs(t), df = df_post) 

    # collect stats
    prot_stats <- data.frame(Contrast=comparison,
			     log2FC=FC, 
			     percentControl=2^FC, 
			     Pvalue=p,
			     Tstatistic=t, 
			     SE=sqrt(variance), 
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
