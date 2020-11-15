#' lmTestContrast
#' @import dplyr
#' @export lmTestContrast 

lmTestContrast <- function(fm, contrast, s2_prior=0, df_prior=0) {

	# check input 
	stopifnot(all(names(contrast) %in% names(coef(fm))))

	# create comparison
	pos_coef <- names(contrast)[contrast>0]
	neg_coef <- names(contrast)[contrast<0]
	comparison <- paste(pos_coef,neg_coef,sep="-")

	# fit the model
	fm <- lm(fx, data = msstats_prot %>% subset(Protein == swip))

	# compute FC
	FC <- as.numeric(LT %*% coeff)
	coeff <- fm$coefficients

	# compute S2 and DF
	av <- anova(fm)
	s2 <- av["Residuals", "Mean Sq"]
	s2_df <- av["Residuals", "Df"]

	# compute posterior S2
	s2_post <- (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)

	# compute variance
	variance <- diag(t(LT) %*% summary(fm)$cov.unscaled %*% LT) * s2_post

	# compute posterior DF
	df_post <- s2_df + df_prior

	# compute t-statistic
	t <- FC / sqrt(variance)

	# compute p-value
	p <- 2 * pt(-abs(t), df = df_post)
	
	# collect stats
	contrast_stats <- data.frame(Contrast=comparison,
					 log2FC=FC, 
					 percentControl=2^FC, 
					 SE=sqrt(variance), 
					 Tstatistic=t, 
					 Pvalue=p,
					 DF=df_post)
	# return the statistics
	return(contrast_stats)
}
