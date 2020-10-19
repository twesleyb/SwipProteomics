#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

# input:
root <- "~/projects/SwipProteomics"
input_prot <- file.path(root,"data","msstats_prot.rda")
input_contrasts <- file.path(root,"data","msstats_contrasts.rda")

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root,quiet=TRUE) }

# imports
suppressPackageStartupMessages({
  library(dplyr)
})


## Functions ------------------------------------------------------------------
# Modified from internal MSstatsTMT functions.

subprot <- function(protein) {
	subdat <- msstats_prot %>% filter(Protein == protein)
	return(subdat)
}

calcPosterior <- function(s2, s2_df, s2.prior, df.prior) {
  # a function to compute s2 posterior
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}


updateLModel <- function(fx, data, mf.final = NULL, ..., change.contr = FALSE) {
  ## PROBLEMATIC
  object <- lmerTest::lmer(fx,data)
  # generate a deviance function as a function of optimal parameters
  #' @importFrom stats formula getCall terms update.formula
  if (is.null(call <- getCall(object))) {
    stop("object should contain a 'call' component")
  }
  # This is where devfun is happening... sneaky arg devFunOnly
  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(mf.final)) {
    call$formula <- update.formula(formula(object), mf.final)
  }
  if (any(grepl("sample", call))) {
    call <- as.list(call)[-which(names(as.list(call)) %in% c("data", "subset"))]
    call[["data"]] <- quote(model.frame(object))
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
    }
  }
  if (change.contr) {
    mm <- model.matrix(object)
    contr <- attr(mm, "contrasts")
    ## change contrasts for F tests calculations
    ## list of contrasts for factors
    if (change.contr && length(which(unlist(contr) != "contr.SAS")) > 0) {
      names.facs <- names(contr)
      l.lmerTest.private.contrast <- as.list(rep("contr.SAS", length(names.facs)))
      names(l.lmerTest.private.contrast) <- names(contr)
      call[["contrasts"]] <- l.lmerTest.private.contrast
    }
    else if (!is.null(contr)) {
      call[["contrasts"]] <- contr[names(contr) %in% attr(terms(call$formula), "term.labels")]
    }
  }
  call <- as.call(call)
  ff <- environment(formula(object))
  pf <- parent.frame() ## save parent frame in case we need it
  sf <- sys.frames()[[1]]
  ff2 <- environment(object)
  tryCatch(eval(call, envir = ff),
    error = function(e) {
      tryCatch(eval(call, envir = sf),
        error = function(e) {
          tryCatch(eval(call, envir = pf),
            error = function(e) {
              eval(call, envir = ff2)
            }
          )
        }
      )
    }
  )
} #EOF


calcApvarcovar <- function(fx,data,thopt,sigma) {
  ## PROBLEMATIC bc calls devFun
  # calc asymptotic variance covariance matrix given params sigma and theta
  fm <- lmerTest::lmer(fx,data)
  dd <- devFun(fx,data)
  h <- hessian(dd, c(thopt, sigma = sigma))
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
  return(A)
} #EOF


devFun <- function(fx,data) {
  ## PROBLEMATIC
  fm <- lmerTest::lmer(fx,data)
  # devfun function as a function of optimal parameters
  #' @importFrom lme4 isGLMM isLMM getME
  #' @importFrom methods is
  # NOTE: appears to be pretty much verbatim from devfun5 in lmertest
  stopifnot(is(fm, "merMod"))
  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }
  n <- nrow(fm@pp$V)
  ff <- updateLModel(fx, data, devFunOnly = TRUE) # calls to updateLModel
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
  class(ans) <- "devFun"
  return(ans)
} #EOF


hessian <- function(fun, x, fx = NULL, delta = 1e-4, ...) {
  # calculate hessian matrix
  # from lmertest::myhess in deriv.R
  nx <- length(x)
  fx <- if (!is.null(fx)) fx else fun(x, ...)
  H <- array(NA, dim = c(nx, nx))
  for (j in 1:nx) {
    ## Diagonal elements:
    xadd <- xsub <- x
    xadd[j] <- x[j] + delta
    xsub[j] <- x[j] - delta
    H[j, j] <- (fun(xadd, ...) - 2 * fx +
      fun(xsub, ...)) / delta^2
    ## Upper triangular (off diagonal) elements:
    for (i in 1:nx) {
      if (i >= j) break
      xaa <- xas <- xsa <- xss <- x
      xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
      xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
      xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
      xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
      H[i, j] <- (fun(xaa, ...) - fun(xas, ...) -
        fun(xsa, ...) + fun(xss, ...)) /
        (4 * delta^2)
    }
  }
  ## Fill in lower triangle:
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  return(H)
} #EOF


vcovLThetaLM <- function(fx,data) {
  ## PROBLEMATIC
  fm <- lmerTest::lmer(fx,data)
  # returns Lc %*% vcov as a function of theta parameters %*% t(Lc)
  #' @importFrom methods is
  #' @importFrom lme4 isGLMM isLMM fixef
  # NOTE: seems to be from vcovJSStheta2 in devfunsLmerTest.R
  stopifnot(is(fm, "merMod"))
  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }
  ff2 <- updateLModel(fx, data, devFunOnly = TRUE)
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
  class(ans) <- "vcovLThetaLM"
  return(ans)
}


gradient <- function(fun, x, delta = 1e-4,
                    method = c("central", "forward", "backward"), ...) {
  # calculate gradient
  # from mygrad in lmertest deriv.R
  method <- match.arg(method)
  nx <- length(x)
  if (method %in% c("central", "forward")) {
    Xadd <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) + diag(delta, nx)
    fadd <- apply(Xadd, 1, fun, ...)
  }
  if (method %in% c("central", "backward")) {
    Xsub <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) - diag(delta, nx)
    fsub <- apply(Xsub, 1, fun, ...) ## eval.parent perhaps?
  }
  res <- switch(method,
    "forward" = (fadd - fun(x, ...)) / delta,
    "backward" = (fun(x, ...) - fsub) / delta,
    "central" = (fadd - fsub) / (2 * delta)
  )
  return(res)
} #EOF


## load the data --------------------------------------------------------------

data(swip)
protein = swip

# load msstats preprocessed protein data from SwipProteomics in root/data
load(file=input_prot) # msstats_prot
load(file=input_contrasts) # msstats_contrasts

# Munge sample annotations
genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
subject <- as.numeric(interaction(msstats_prot$Mixture,msstats_prot$Genotype))
msstats_prot$Genotype <- factor(genotype,levels=c("Control","Mutant"))
msstats_prot$BioFraction <- factor(biofraction,levels=c("F4","F5","F6","F7","F8","F9","F10"))
msstats_prot$Subject <- as.factor(subject)

fx <- formula("Abundance ~ 1 + (1|Mixture) + BioFraction + Genotype")
fm <- lmerTest::lmer(fx,data=subprot(swip))

library(merDeriv)

#s2 = sigma_squared from rho
av <- anova(fm)
s2 <- av$"Mean Sq" / av$"F value"

# Beta = estimates(lmer,Condition)
# variance-covariance matrix of fixed effects (Condition)
#V_hat = merDeriv::vcov.lmerMod(fm,full=FALSE) 
