#' lmerTestProtein
#' @import lme4 lmerTest
#' @export lmerTestProtein

lmerTestProtein <- function(protein, fx, msstats_prot, contrasts, gof=FALSE) {
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
				  nm=c("R2m_fixef", "R2c_total"))
	  r2_nakagawa["R2r_mixef"] <- r2_nakagawa[2] - r2_nakagawa[1]
	  r2_nakagawa <- r2_nakagawa[c(3,1,2)]
  }
  # compute Satterthwaite degrees of freedom and other key statistics
  model_summary <- summary(fm, ddf = "Satterthwaite")
  s2_df <- as.numeric(model_summary$coefficients[,"df"][1]) 
  coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)
  sigma <- model_summary$sigma # == stats::sigma(fm)
  theta <- model_summary$optinfo$val # aka thopt == lme4::getME(fm, "theta")
  vcov <- model_summary$vcov # variance-covariance matrix
  # calcuate symtoptic var-covar matrix
  A <- fm@vcov_varpar
  ## Loop through contrasts to perform tests.
  ## FIXME: should be a function that does one iter
  stats_list <- list()
  for (contrast in contrasts) {
    # we compute the unmoderated statistics:
    df_prior=0
    s2_prior=0
    # things that depend upon the contrast
    se2 <- as.numeric(contrast %*% vcov %*% contrast) # == variance 
    # calculate posterior s2
    s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)
    # calculate gradient
    g <- c(contrast %*% fm@Jac_list[[1]] %*% contrast,
	   contrast %*% fm@Jac_list[[2]] %*% contrast)
    denom <- as.numeric(t(g) %*% A %*% g)
    # NOTE: which is correct?
    #df_post <- (2 * se2) / (denom + df_prior) # or se2^2 ???
    df_post <- (2 * se2^2) / (denom + df_prior)
    # compute fold change and the t-statistic
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(se2) 
    # compute the p-value given t-statistic and df.post
    p <- 2 * pt(-abs(t), df = df_post) 
    comparison <- paste(names(contrast)[contrast == +1], 
		        names(contrast)[contrast == -1],sep="-")
    stats_list[[comparison]] <- data.frame(protein=protein,contrast=comparison,
  		   log2FC=FC, percentControl=2^FC, Pvalue=p,
  		   Tstatistic=t, SE=sqrt(se2), DF=df_post, 
		   isSingular=lme4::isSingular(fm))
  } # EOL through contrasts
  ## compile results
  rho <- list()
  rho$protein <- protein
  rho$formula <- fx
  rho$stats <- bind_rows(stats_list)
  if (gof) {
	  rho$gof <- r2_nakagawa
  }
  return(rho)
} #EOF
