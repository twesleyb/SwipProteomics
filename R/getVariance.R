#' export getVariance

getVariance <- function(fm) {
  # Get the fixed and mixed effects variance componets of the model
  var_df <- as.data.frame(lme4::VarCorr(fm))
  var_mixef <- setNames(var_df$vcov, nm = var_df$grp)
  var_fixef <- stats::var(as.vector(lme4::fixef(fm) %*%
    t(lme4::getME(fm, "X"))))
  variance <- c(var_mixef, "Fixed" = var_fixef)
  return(variance)
}
