#!/usr/bin/env Rscript

.vcovLThetaL <- function(fm) {
#############################################
## returns Lc %*% vcov as a function of theta parameters %*% t(Lc)
#############################################
#' @importFrom lme4 isGLMM isLMM fixef
#' @importFrom methods is
#' @keywords internal
	# NOTE: seems to be from vcovJSStheta2 in devfunsLmerTest.R
  stopifnot(is(fm, "merMod"))

  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }

  ff2 <- .updateModel(fm, devFunOnly = TRUE)

  envff2 <- environment(ff2)

  if (lme4::isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)

      sigma2 <- thpars[np]^2
      ff2(thpars[-np])

      vcov_unscaled <- tcrossprod(envff2$pp$RXi())
      vcov_out <- sigma2 * vcov_unscaled

      return(list(
        varcor = as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)),
        unscaled.varcor = vcov_unscaled,
        sigma2 = sigma2
      ))
    }
  }
  class(ans) <- ".vcovLThetaL"

  ans
}
