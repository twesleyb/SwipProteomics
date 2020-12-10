#' getVariance
#'
#' extract the fixed and mixed variance components from lmer
#'
#' @export getVariance
#'
#' @importFrom lme4 VarCorr fixef getME

getVariance <- function(fm) {
  # Get the fixed and mixed effects variance components of the model

  var_df <- as.data.frame(lme4::VarCorr(fm))
  var_mixef <- setNames(var_df$vcov, nm = var_df$grp)

  coeff <- lme4::fixef(fm)
  var_fixef <- stats::var(as.vector(coeff %*% t(lme4::getME(fm, "X"))))

  variance <- c(var_mixef, "Fixed" = var_fixef)

  return(variance)
}
