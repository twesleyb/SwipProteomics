#!/usr/bin/env Rscript

.calcApvar <- function(rho) {
#############################################
## calculates asymptotic variance covariance matrix of variance parameters based on theta
#############################################
#' @keywords internal
  ## based on theta parameters and sigma
  # NOTE: not from lmertest
  dd <- .devfunTheta(rho$model)
  h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma))

  ch <- try(chol(h), silent = TRUE)
  if (inherits(ch, "try-error")) {
    return(rho)
  }
  A <- 2 * chol2inv(ch)

  eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  isposA <- TRUE
  if (min(eigval) < sqrt(.Machine$double.eps)) { ## tol ~ sqrt(.Machine$double.eps)
    isposA <- FALSE
  }

  if (!isposA) {
    print("Asymptotic covariance matrix A is not positive!")
  }

  A
}
