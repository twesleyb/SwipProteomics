#!/usr/bin/env Rscript

#' calcApvar
#' calculate the asymptotic variance-covariance matrix of variance parameters 
#' @keywords internal

.calcApvar <- function(fm,thopt,sigma) {
  # based on theta parameters and sigma

  dd <- .devfunTheta(fm)
  h <- .myhess(dd, c(thopt, sigma))
  ch <- try(chol(h), silent = TRUE)

  if (inherits(ch, "try-error")) {
    return(rho)
  }

  A <- 2 * chol2inv(ch)

  eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values

  if (min(eigval) < sqrt(.Machine$double.eps)) { # tolerance
    warning("Asymptotic covariance matrix A is not positive!")
  }

  A
}
