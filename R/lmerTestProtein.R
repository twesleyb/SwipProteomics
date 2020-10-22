#' lmerTestProtein
#' @importFrom lmerTest
#' @export lmerTestContrast lmerTestProtein

lmerTestContrast <- function(contrast, fm, df_prior=0,s2_prior=0) {

    # define the comparison to be tested
    pos_group <- names(contrast[contrast==1])
    neg_group <- names(contrast[contrast==-1])
    comparison <- paste(pos_group,neg_group,sep="-")

    # get model's summary and compute df
    model_summary <- summary(fm, ddf = "Satterthwaite")

    # get df cooresponding to coeff in comparison -- but which coeff?
    s2_df <- model_summary$coefficients[,"df"][c(pos_group,neg_group)][1]
    coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)

    # extract sigma, theta, and variance-covariance matrix
    sigma <- model_summary$sigma # == stats::sigma(fm)
    theta <- model_summary$optinfo$val # aka thopt == lme4::getME(fm, "theta")
    vcov <- model_summary$vcov # variance-covariance matrix

    # extract asymtoptic var-covar matrix from fit model
    A <- fm@vcov_varpar

    # vcov_varpara numeric matrix holding the asymptotic variance-covariance 
    # matrix of the variance parameters (including sigma).
    # we compute the unmoderated statistics:

    # things that depend upon the contrast
    se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance 

    # calculate posterior s2
    s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)

    # calculate gradient
    g <- c(contrast %*% fm@Jac_list[[1]] %*% contrast,
	   contrast %*% fm@Jac_list[[2]] %*% contrast)

    # Jac_list - a list of gradient matrices (Jacobians) for the gradient of the
    # variance-covariance of beta with respect to the variance parameters, where
    # beta are the mean-value parameters available in fixef(object).

    # given gradient and asymptoptic var-covar, calculate denom (?)
    denom <- as.numeric(t(g) %*% A %*% g)

    # compute posterior df
    # bayesian? posterior proporional to inverse of - σ² ?
    #df_post <- (2 * se2) / (denom + df_prior)
    df_post <- se2 / (denom + df_prior)

    # compute fold change and the t-statistic
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(se2) # 

    # compute the p-value given t-statistic and df.post
    p <- 2 * pt(-abs(t), df = df_post) 

    # collect stats
    stats_df <- data.frame(contrast=comparison,
			   log2FC=FC, 
			   percentControl=2^FC, 
			   Pvalue=p,
			   Tstatistic=t, 
			   SE=sqrt(se2), 
			   DF=df_post, 
			   isSingular=lme4::isSingular(fm))

    return(stats_df)
} #EOF
#' lmerTestProtein
#' @import lme4 lmerTest
#' @export lmerTestProtein

lmerTestProtein <- function(protein, fx, msstats_prot, contrasts, gof=FALSE) {

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

  # evaluate R2 for marginal (fixed) and conditional (total) effects
  if (gof) {
	  r2_nakagawa <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),
				  nm=c("R2_fixef", "R2_total"))
	  r2_nakagawa["R2_mixef"] <- r2_nakagawa[2] - r2_nakagawa[1]
	  r2_nakagawa <- r2_nakagawa[c(3,1,2)]
  }

  # evaluate statistical comparisons between contrasts
  stats <- lapply(contrasts, lmerTestContrast, fm)

  # compile results
  rho <- list()
  rho$protein <- protein
  rho$formula <- fx
  rho$stats <- as.data.frame(dplyr::bind_rows(stats))
  if (gof) {
	  rho$gof <- r2_nakagawa
  }

  return(rho)
} #EOF
