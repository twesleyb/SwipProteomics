#' lmerTestProtein
#' @import lme4 lmerTest
#' @export lmerTestProtein

lmerTestProtein <- function(protein, fx, msstats_prot, 
			    contrast_matrix, s2_prior=0, df_prior=0) {
  # subset the data
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
	# the data cannot contain missing values
  	warning("The data cannot contain missing values.")
        return(NULL)
  }
  # fit the model
  fm <- lmerTest::lmer(fx, data=subdat)
  # compute Satterthwaite degrees of freedom and other key statistics
  model_summary <- summary(fm, ddf = "Satterthwaite")
  s2_df <- as.numeric(model_summary$coefficients[,"df"][1]) 
  coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)
  sigma <- model_summary$sigma # == stats::sigma(fm)
  theta <- model_summary$optinfo$val # aka thopt == lme4::getME(fm, "theta")
  vcov <- model_summary$vcov # variance-covariance matrix
  se2 <- as.numeric(contrast_matrix %*% vcov %*% contrast_matrix) # == variance
  s2_prior = sigma^2
  # calculate posterior s2
  s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)
  # calcuate symtoptic var-covar matrix
  A <- fm@vcov_varpar
  # FIXME: there might not always be two?
  g <- c(contrast_matrix %*% fm@Jac_list[[1]] %*% contrast_matrix,
	 contrast_matrix %*% fm@Jac_list[[2]] %*% contrast_matrix)
  denom <- as.numeric(t(g) %*% A %*% g)
  # NOTE: which is correct?
  df_post <- (2 * se2) / (denom + df_prior) # or se2^2 ???
  #df_post <- (2 * se2^2) / (denom + df_prior)
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% coeff)[, 1]
  t <- FC / sqrt(se2) 
  # compute the p-value given t-statistic and df.post
  p <- 2 * pt(-abs(t), df = df_post) 
  # compile results
  rho <- list()
  rho$protein <- protein
  rho$model <- fx
  comparison <- paste(names(contrast_matrix)[contrast_matrix == +1], 
		      names(contrast_matrix)[contrast_matrix == -1],sep="-")
  rho$stats <- data.frame(protein=protein,contrast=comparison,
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
  		 Tstatistic=t, SE=sqrt(se2), DF=df_post, 
		 isSingular=lme4::isSingular(fm))
  return(rho)
} #EOF
