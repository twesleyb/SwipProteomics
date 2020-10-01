#!/usr/bin/env Rscript

.devfunTheta <- function(fm) {
#############################################
## devfun function as a function of optimal parameters
#############################################
#' @importFrom lme4 isGLMM isLMM getME
#' @importFrom methods is
#' @keywords internal
  # NOTE: appears to be pretty much verbatim from devfun5 in lmertest
  stopifnot(is(fm, "merMod"))

  np <- length(fm@pp$theta)
  nf <- length(fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }
  n <- nrow(fm@pp$V)

  ff <- .updateModel(fm, devFunOnly = TRUE)
  reml <- lme4::getME(fm, "is_REML")

  envff <- environment(ff)

  if (lme4::isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)

      message("What is thpars?")
      print(class(thpars))
      print(thpars)

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
