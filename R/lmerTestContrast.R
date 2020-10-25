#' lmerTestContrast 
#' @import dplyr lmerTest data.table
#' @export lmerTestContrast


lmerTestContrast <- function(fm, contrast, 
			     df_prior=0,s2_prior=0) {

    # comparison should be a numeric and sum of vector == 0
    stopifnot(inherits(contrast,"numeric"))
    stopifnot(sum(contrast[contrast<0],contrast[contrast>0])==0)
    pos_group <- names(contrast[contrast>0])
    neg_group <- names(contrast[contrast<0])
    comparison <- paste(pos_group,neg_group,sep="-")

    # compute Satterthwaite degrees of freedom
    model_summary <- summary(fm, ddf = "Satterthwaite")

    # collect model's df and coefficients
    s2_df <- as.numeric(model_summary$coefficients[,"df"][pos_group])[1]
    coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)

    # extract sigma and theta
    sigma <- as.numeric(model_summary$sigma) # == stats::sigma(fm)
    theta <- as.numeric(model_summary$optinfo$val) # == lme4::getME(fm, "theta")

    # calculate posterior s2
    s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)

    # compute the unscaled covar matrix with unsc()
    unscaled_vcov <- fm@pp$unsc()

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

    # compute fold change and the t-statistic
    FC <- (contrast %*% coeff)[, 1]
    t <- FC / sqrt(se2) 

    # compute the p-value given t-statistic and posterior degrees of freedom
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
