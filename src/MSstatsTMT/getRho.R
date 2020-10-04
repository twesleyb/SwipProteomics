#!/usr/bin/env Rscript

#' @importFrom lme4 fixef getME
#' @importFrom stats sigma
#' @keywords internal

#' @export

getRho <- function(fm) {
  # return a list rho with statistics from the fit model
  # added by TWAB
  rho <- list()

  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$thopt <- lme4::getME(fm, "theta")

  # calculate degrees of freedom and sigma^2:
  av <- anova(fm)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"

  # calcuate symtoptic var-covar matrix
  rho$model <- fm
  rho$A <-.calcApvar(rho)

  return(rho)
}
