#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: 

# input:
# * msstats_prot.rda

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
})


## MSstatsTMT internal functions ----------------------------------------------

calcPosterior <- function(s2, s2_df, s2.prior = 0, df.prior = 0) {
  # A function to compute s2 posterior
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}


updateLModel <- function(object, mf.final = NULL, ..., change.contr = FALSE) {
  ## devfun function as a function of optimal parameters
  #' @importFrom stats formula getCall terms update.formula
  #' @keywords internal
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


calcApvarcovar <- function(rho) {
  # calc asymptotic variance covariance matrix of variance parameters sigma 
  # and theta
  dd <- devFun(rho$model)
  h <- myhessian(dd, c(rho$thopt, sigma = rho$sigma))
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


devFun <- function(fm) {
  ## devfun function as a function of optimal parameters
  #' @importFrom lme4 isGLMM isLMM getME
  #' @importFrom methods is
  #' @keywords internal
  # NOTE: appears to be pretty much verbatim from devfun5 in lmertest
  stopifnot(is(fm, "merMod"))
  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) {
    np <- np + 1L
  }
  n <- nrow(fm@pp$V)
  ff <- updateLModel(fm, devFunOnly = TRUE)
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


myhessian <- function(fun, x, fx = NULL, delta = 1e-4, ...) {
  ## calculate hessian matrix
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


vcovLThetaLM <- function(fm) {
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
  ff2 <- updateLModel(fm, devFunOnly = TRUE)
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


mygradient <- function(fun, x, delta = 1e-4,
                    method = c("central", "forward", "backward"), ...) {
  ## calculate gradient
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

# load msstats preprocessed protein data from SwipProteomics in root/data
devtools::load_all()
data(msstats_prot) 

# Munge sample annotations
genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
msstats_prot$Genotype <- as.factor(genotype)
msstats_prot$BioFraction <- biofraction


## Analysis of Washc4/Swip -----------------------------------------------------

# Swip/Washc4
swip = protein = "Q3UMB9"

# NOTE: the moderated analysis requires all models be fit. Here we demonstrate
# fitting just a single protein.
# moderated = FALSE

# subset the data for a given protein [START]
subdat <- msstats_prot %>% filter(Protein == protein)

# The linear mixed effects model to be fit:
fx1 <- formula("Abundance ~ 0 + (1|BioFraction) + Genotype")

# fit the model:
fm1 <- lmerTest::lmer(fx1, data=subdat)
message("fit lmer: ", fx1)

# examine the fit:
print(summary(fm1))

# compute some model statistics, store in list rho
rho <- list()
rho$protein <- protein
rho$formula <- fx1
rho$model <- fm1
rho$data <- subdat

# here we compute the unmoderated statistics
rho$df.prior <- 0
rho$s2.prior <- 0

# check if singular (boundary) fit
rho$isSingular <- lme4::isSingular(fm1)

# calculate coeff, sigma, and theta:
rho$coeff <- lme4::fixef(fm1)
rho$sigma <- stats::sigma(fm1)
rho$thopt <- lme4::getME(fm1, "theta")

# calculate degrees of freedom and sigma^2:
av <- anova(fm1)
rho$s2_df <- av$DenDF
rho$s2 <- av$"Mean Sq" / av$"F value"

# calcuate symtoptic var-covar matrix
# NOTE: utlilizes MSstatsTMT internal function .calcApvar to do calculation
rho$A <- calcApvarcovar(rho) 

# define contrast matrix for comparison to control
cm1 <- rho$coeff
cm1[1] <- -1
cm1[2] <- 1
rho$contrasts <- cm1

# we store the proteins statistics in a list
stats_list <- list()

# calculate posterior s2
s2.post <- calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)

# compuate variance-covariance matrix
# NOTE: uses vcovLThetaLM
vss <- vcovLThetaLM(fm1)
varcor <- vss(t(cm1), c(rho$thopt, rho$sigma))
vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
se2 <- as.numeric(t(cm1) %*% as.matrix(vcov) %*% cm1)

# calculate variance
vcov.post <- varcor$unscaled.varcor * s2.post
variance <- as.numeric(t(cm1) %*% as.matrix(vcov.post) %*% cm1)

# calculate degrees of freedom
# given params theta and sigma from lmer and the Acovar
# NOTE: uses mygradient
g <- mygradient(function(x) vss(t(cm1), x)$varcor, c(rho$thopt, rho$sigma))
denom <- t(g) %*% rho$A %*% g

# compute df.posterior
df.post <- 2 * (se2)^2 / denom + rho$df.prior

# compute fold change and the t-statistic
FC <- (cm1 %*% rho$coeff)[, 1]
t <- FC / sqrt(variance)

# compute the p-value
p <- 2 * pt(-abs(t), df = df.post)

# munge comparison
comparison <- gsub("Genotype","",paste(rev(names(cm1)),collapse="-"))

# compile results
df <- data.frame(protein=protein,
		 contrast=comparison,
		 log2FC=FC,
		 percentControl=2^FC,
		 Pvalue=p,
		 Tstatistic=t,
		 SE=sqrt(variance),
		 DF=df.post,
		 isSingular=rho$isSingular)

knitr::kable(df)


## Declare a function that does all the work of above -------------------------

# * msstats_prot
# * protein
# * fm1 declared above (lmer formula)
# * cm1 declared above (contrast matrix)

lmerTestProtein <- function(msstats_prot,protein,fm1,cm1) {
  # perfrom the lmer-based testing of conditioned means for a given protein
  # subset the data for a given protein
  subdat <- msstats_prot %>% filter(Protein == protein)
  # the data cannot contain missing values
  if (any(is.na(subdat))) {
  	stop("The data cannot contain missing values.")
  }
  # fit the model:
  fm1 <- suppressMessages({try(lmerTest::lmer(fx1, data=subdat),silent=TRUE)})
  # FIXME: catch errors/warnings!
  # compute some model statistics, store in list rho
  rho <- list()
  rho$protein <- protein
  rho$formula <- fx1
  rho$model <- fm1
  rho$data <- subdat
  # here we compute the unmoderated statistics
  rho$df.prior <- 0
  rho$s2.prior <- 0
  # check if singular (boundary) fit
  rho$isSingular <- lme4::isSingular(fm1)
  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm1)
  rho$sigma <- stats::sigma(fm1)
  rho$thopt <- lme4::getME(fm1, "theta")
  # calculate degrees of freedom and sigma^2:
  av <- anova(fm1)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"
  # calcuate symtoptic var-covar matrix
  # NOTE: utlilizes MSstatsTMT internal function .calcApvar to do calculation
  rho$A <- calcApvarcovar(rho) 
  # we store the proteins statistics in a list
  stats_list <- list()
  # calculate posterior s2
  s2.post <- calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)
  # compuate variance-covariance matrix
  vss <- vcovLThetaLM(fm1)
  varcor <- vss(t(cm1), c(rho$thopt, rho$sigma))
  vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
  se2 <- as.numeric(t(cm1) %*% as.matrix(vcov) %*% cm1)
  # calculate variance
  vcov.post <- varcor$unscaled.varcor * s2.post
  variance <- as.numeric(t(cm1) %*% as.matrix(vcov.post) %*% cm1)
  # calculate degrees of freedom
  # given params theta and sigma from lmer and the Acovar
  # NOTE: careful formatting can break things
  g <- mygradient(function(x)vss(t(cm1), x)$varcor, c(rho$thopt, rho$sigma))
  denom <- t(g) %*% rho$A %*% g
  # compute df.posterior
  df.post <- 2 * (se2)^2 / denom + rho$df.prior
  # compute fold change and the t-statistic
  FC <- (cm1 %*% rho$coeff)[, 1]
  t <- FC / sqrt(variance)
  # compute the p-value
  p <- 2 * pt(-abs(t), df = df.post)
  # munge comparison
  comparison <- gsub("Genotype","",paste(rev(names(cm1)),collapse="-"))
  # compile results
  df <- data.frame(protein=protein,
  		 contrast=comparison,
  		 log2FC=FC,
  		 percentControl=2^FC,
  		 Pvalue=p,
  		 Tstatistic=t,
  		 SE=sqrt(variance),
  		 DF=df.post)
  return(df)
} #EOF

#------------------------------------------------------------------------------
## BLARGH FAILS: first protein error in t(g) rho$A must be numeric, but runs
## above probably because of calcApvar within a loop

#            |-- devfunTheta ... fm@ updateModel (weird env stuff)
# calcApvar -|-- myhess ... lmertest::mytest() (deriv.R)

# check:
result <- lmerTestProtein(msstats_prot,swip,fm1,cm1)

# Loop through all proteins
results_list <- list()
proteins <- unique(as.character(msstats_prot$Protein))
message("\nAnalyzing all proteins.")

# ERROR:
for (i in seq(proteins)) {
  if (i==1) {
	# init pbar
	pbar <- txtProgressBar(max=length(proteins),style=3)
  }
  results_list[[i]] <- lmerTestProtein(msstats_prot,proteins[i],fm1,cm1)
  setTxtProgressBar(pbar,match(protein,proteins))
}
close(pbar)
