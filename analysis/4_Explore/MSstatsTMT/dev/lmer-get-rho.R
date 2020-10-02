#!/usr/bin/env Rscript

#' @importFrom lme4 fixef getME
#' @importFrom stats sigma
#' @keywords internal

get_rho <- function(fm) {
  rho <- list()
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$theta <- lme4::getME(fm, "theta")

  av <- anova(fm)
  rho$df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"

  # add apvar
  A <-.calcApvar(rho)

  return(rho)
}
