#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

## prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root); devtools::load_all(root)

## inputs:
data(msstats_prot)
#data(gene_map)
data(swip)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
})


## Functions ------------------------------------------------------------------
# Modified from internal MSstatsTMT functions.

# * calcPosterior
# * updateLModel
# * calcApvarcovar
# * devFun
# * vcovLThetaLM
# * myhessian
# * mygradient
# * lmerTestProtein - model fitting and statistical testing for a given protein

# NOTE: I added fx, data to all functions that handle the lmer model. 
# These functions init a fit model by fitting the function to the data.
# NOTE: Names have been changed by MSstats and then again by myself for clarify
# or uniqueness.


calcPosterior <- function(s2, s2_df, s2.prior, df.prior) {
  # a function to compute s2 posterior
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}


updateLModel <- function(fx, data, mf.final = NULL, ..., change.contr = FALSE) {
  object <- lmerTest::lmer(fx,data)
  # NOTE: from lmerTest update
  # generate a deviance function as a function
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


#calcApvarcovar <- function(fx,data,thopt,sigma) {
#  # why not:
#  # as_lmerModLmerTest(fx, tol = 1e-08)
#  # vcov_varpar of param theta sigma
#  # calc asymptotic variance covariance matrix given params sigma and theta
#  #' @importFrom lmerTest lmer
#  fm <- lmerTest::lmer(fx,data)
#  dd <- devFun(fx,data)
#  h <- myhessian(dd, c(thopt, sigma = sigma))
#  ch <- try(chol(h), silent = TRUE)
#  if (inherits(ch, "try-error")) {
#    return(rho)
#  }
#  A <- 2 * chol2inv(ch)
#  eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
#  isposA <- TRUE
#  if (min(eigval) < sqrt(.Machine$double.eps)) { ## tol ~ sqrt(.Machine$double.eps)
#    isposA <- FALSE
#  }
#  if (!isposA) {
#    print("Asymptotic covariance matrix A is not positive!")
#  }
#  return(A)
#} #EOF


devFun <- function(fx,data) {
  # NOTE: requires updateLModel
  fm <- lmerTest::lmer(fx,data)
  # devfun function as a function of optimal parameters
  #' @importFrom methods is
  #' @importFrom lmerTest lmer
  #' @importFrom lme4 isGLMM isLMM getME
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


#myhessian <- function(fun, x, fx = NULL, delta = 1e-4, ...) {
#  # calculate hessian matrix
#  # from lmertest::myhess in deriv.R
#  nx <- length(x)
#  fx <- if (!is.null(fx)) fx else fun(x, ...)
#  H <- array(NA, dim = c(nx, nx))
#  for (j in 1:nx) {
#    ## Diagonal elements:
#    xadd <- xsub <- x
#    xadd[j] <- x[j] + delta
#    xsub[j] <- x[j] - delta
#    H[j, j] <- (fun(xadd, ...) - 2 * fx +
#      fun(xsub, ...)) / delta^2
#    ## Upper triangular (off diagonal) elements:
#    for (i in 1:nx) {
#      if (i >= j) break
#      xaa <- xas <- xsa <- xss <- x
#      xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
#      xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
#      xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
#      xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
#      H[i, j] <- (fun(xaa, ...) - fun(xas, ...) -
#        fun(xsa, ...) + fun(xss, ...)) /
#        (4 * delta^2)
#    }
#  }
#  ## Fill in lower triangle:
#  H[lower.tri(H)] <- t(H)[lower.tri(H)]
#  return(H)
#} #EOF


vcovLThetaLM <- function(fx,data) {
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
      # NOTE: RXi is a function that returns the inverse of the dense 
      # downdated Cholesky factor for the fixed-effects parameters
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


lmerTestProtein <- function(protein, fx, msstats_prot, contrast_matrix) {
  # slimed down, but requires:
  # * calcPosterior
  # * vcovLThetaLM(fx,subdat)
  # * mygradient(myfun, c(theta, sigma))
  # here we compute the unmoderated statistics
  df_prior <- 0
  s2_prior <- 0
  # subset the data
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
	# the data cannot contain missing values
  	warning("The data cannot contain missing values.")
        return(NULL)
  }
  # fit the model
  fm <- lmerTest::lmer(fx, data=subdat)
  # compute Satterthwaite degrees of freedom and other key statistics
  model_summary <- summary(fm, ddf = "Satterthwaite")
  s2_df <- as.numeric(model_summary$coefficients[,"df"][1])
  coeff <- model_summary$coeff[,"Estimate"]
  sigma <- model_summary$sigma # == sigma(fm)
  theta <- model_summary$optinfo$val # aka thopt
  vcov <- model_summary$vcov # variance-covariance matrix
  se2 <- as.numeric(contrast_matrix %*% vcov %*% contrast_matrix) # variance
  # calculate posterior s2
  s2_post <- calcPosterior(sigma^2, s2_df, s2_prior, df_prior)
  # calculate degrees of freedom posterior
  # given params theta and sigma from lmer and the covariance matrix
  # calcuate symtoptic var-covar matrix
  vss <- vcovLThetaLM(fx,subdat)
  varcor <- vss(t(contrast_matrix), c(theta, sigma))
  vcov.post <- varcor$unscaled.varcor * s2_post
  A <- fm@vcov_varpar
  myfun <- function(x) { vss(t(contrast_matrix), x)$varcor }
  g <- mygradient(myfun, c(theta, sigma))
  denom <- as.numeric(t(g) %*% A %*% g)
  #df.post <- 2 * se2^2 / denom + df.prior # is se2 ^2 correct?
  df_post <- 2 * se2 / denom + df_prior
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% coeff)[, 1]
  t <- FC / sqrt(se2) 
  # compute the p-value given t-statistic and df.post
  p <- 2 * pt(-abs(t), df = df_post) 
  # compile results
  rho <- list()
  rho$protein <- protein
  rho$model <- fx
  comparison <- paste(names(contrast_matrix)[contrast_matrix == +1], 
		      names(contrast_matrix)[contrast_matrix == -1],sep="-")
  rho$stats <- data.frame(protein=protein,contrast=comparison,
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
  		 Tstatistic=t, SE=sqrt(se2), DF=df_post, 
		 isSingular=lme4::isSingular(fm))
  return(rho)
} #EOF


## check Swip's fit -----------------------------------------------------------

# demonstrate fit:
fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

model_summary <- summary(fm0, ddf = "Satterthwaite")
model_summary

# FIXME: need to replace with more reproducible code
knitr::kable(r.squaredGLMM.merMod(fm0))

# build a contrast matrix:
cm0 <- lme4::fixef(fm0)
cm0[] <- 0
cm0["ConditionControl.F7"] <- -1
cm0["ConditionMutant.F7"] <- +1 

# test a comparison defined by contrast_matrix
model0 <- lmerTestProtein(swip, fx0, msstats_prot, cm0)

model0$stats %>% knitr::kable()

## repeat for comparisons across all fractions:

# fit a model
fx1 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Subject)")
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein == swip))
summary(fm1)

knitr::kable(r.squaredGLMM.merMod(fm1))

# generate contrast
cm1 <- lme4::fixef(fm1)
cm1[] <- 0
cm1["GenotypeControl"] <- -1
cm1["GenotypeMutant"] <- +1

# test the model
model1 <- lmerTestProtein(swip, fx1, msstats_prot, cm1)
model1$stats %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

#n_cores <- parallel::detectCores()
#BiocParallel::register(BiocParallel::SnowParam(n_cores))
#
#prots = unique(as.character(msstats_prot$Protein))
#
#results_list <- foreach(protein = prots) %dopar% {
#	suppressMessages({
#	  try(lmerTestProtein(protein, fx1, msstats_prot, cm1),silent=TRUE)
#	})
#} # EOL
#
#
### collect results ------------------------------------------------------------
## collect stats
#idx <- unlist(sapply(results_list,class)) != "try-error"
#filt_list <- results_list[which(idx)]
#results_df <- bind_rows(sapply(filt_list,"[","stats"))
#
## drop singular
#results_df <- results_df %>% filter(!isSingular)
#results_df$isSingular <- NULL
#
## annotate with gene symbols
#idx <- match(results_df$protein,gene_map$uniprot)
#results_df <- tibble::add_column(results_df,
#  				 symbol=gene_map$symbol[idx],
#  				 .after="protein")
#
## adjust pvals 
#results_df <- tibble::add_column(results_df, 
#			 Padjust=p.adjust(results_df$Pvalue,"BH"),
#			 .after="Pvalue")
#
## sort
#results_df <- results_df %>% arrange(Pvalue)
#
## examine top results
#results_df %>% head() %>% knitr::kable()
