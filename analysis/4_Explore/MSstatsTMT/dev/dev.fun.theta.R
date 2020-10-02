#!/usr/bin/env Rscript

#' @importFrom lme4 isGLMM isLMM getME
#' @importFrom methods is
#' @keywords internal

.devfunTheta <- function(fm) {
  # devfun function as a function of optimal parameters
  # NOTE: mostly from devfun5 in lmertest, substitute .updateModel
  # NOTE: try removing .updateModel dependency

  stopifnot(is(fm, "merMod"))

  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }
  n <- nrow(fm@pp$V)

  # here devFunOnly is passed as sneaky arg to .updateModel
  # seems like cutting out updateModel with update doesn't affect things
  #ff <- .updateModel(fm, devFunOnly = TRUE) # does this substituion matter?
  ff <- update(fm, devFunOnly=TRUE)
  reml <- lme4::getME(fm, "is_REML")

  envff <- environment(ff)

  if (lme4::isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)

      ff(thpars[-np])

      sigsq <- thpars[np]^2

      dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1)) / sigsq + n * log(2 * pi * sigsq)

      if (reml) {
        p <- ncol(envff$pp$RX())
        dev <- dev + 2 * determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)
      }
      return(dev)
    }
  }

  attr(ans, "thopt") <- fm@pp$theta
  class(ans) <- ".devfunTheta"
  ans
}
