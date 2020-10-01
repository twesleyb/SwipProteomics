#!/usr/bin/env Rscript

.rhoInit <- function(rho, model, change.contr = FALSE, mf.final = NULL) {
##########################################################################
## Create rho vector containing info about mixed model ##################
##########################################################################
#' @importFrom lme4 fixef getME
#' @importFrom stats sigma
#' @keywords internal
#   https://github.com/rforge/lmertest/blob/master/pkg/lmerTest/R
  if (change.contr) {
    rho$model <- .updateModel(model, mf.final = mf.final, change.contr = change.contr)
  } else {
    rho$model <- model
  }
  rho$fixEffs <- lme4::fixef(rho$model)
  rho$sigma <- stats::sigma(rho$model)
  rho$thopt <- lme4::getME(rho$model, "theta")
  return(rho)
}

#get_rho <- function(fit) {
#	rho <- list()
#	rho[["fixEffs"]] <- lme4::fixef(fit)
#        rho[["sigma"]] <- stats::sigma(rho$model)
#        rho[["thopt"]] <- lme4::getME(rho$model, "theta")
#        return(rho)
#}
