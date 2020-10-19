#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

# input:
root <- "~/projects/SwipProteomics"
# * msstats_prot

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root,quiet=TRUE) }

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
})


## Functions ------------------------------------------------------------------
# Modified from internal MSstats functions.

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
  #' @importFrom lmerTest lmer
  fm <- lmerTest::lmer(fx,data)
  dd <- devFun(fx,data)
  h <- myhessian(dd, c(thopt, sigma = sigma))
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


myhessian <- function(fun, x, fx = NULL, delta = 1e-4, ...) {
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


## load the data --------------------------------------------------------------

# load msstats preprocessed protein data from SwipProteomics in root/data
#devtools::load_all()
load(file.path(root,"data","msstats_prot.rda"))
load(file.path(root,"data","swip.rda")); protein = swip
data(gene_map)

# Munge sample annotations:
# * create Genotype column
# * create BioFraction column
genotype <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",1)
biofraction <- sapply(strsplit(as.character(msstats_prot$Condition),"\\."),"[",2)
subject <- interaction(msstats_prot$Mixture,genotype)
msstats_prot$Genotype <- as.factor(genotype)
msstats_prot$BioFraction <- biofraction
msstats_prot$Subject <- subject

# define protein to be analyzed
proteins = unique(as.character(msstats_prot$Protein))

lmerTestProtein <- function(protein) {
  #############################################################################
  ## Q1. Which model is the correct model?
  ## We are interested in differences between levels of Genotype
  ## adjusted for 'BioFraction'.
  ## [1] all covariates
  # boundary fit -- problem with combination of mixture and subject
  #fx <- formula("Abundance ~ BioFraction + (1|Mixture) + (1|Subject) + Genotype") 
  ## [2] Only mixed effect of Subject:
  fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Subject)")
  # NOTE: what is meaning of intercept term 0 or 1 or just Abundance ~
  # BioFraction... 
  #r.squaredGLMM.merMod(fm)
  #    R2m       R2c
  # [1,] 0.9162647 0.9663459 <-- better fit(?) (more var explained)
  # R2c is total variance explained
  # R2m is the variance explained by fixed effects
  # remaineder = ~0.05 is the variance explained by mixed effects.
  ## [3] Only mixed effect of Mixture:
  #fx <- formula("Abundance ~ BioFraction + (1|Mixture) + Genotype")
  #r.squaredGLMM.merMod(fm)
  #    R2m       R2c
  #[1,] 0.9222618 0.9355818  <-- R2c = total variance explained 
  #############################################################################
  # the data cannot contain missing values
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
  	warning("The data cannot contain missing values.")
        return(NULL)
  }
  fm <- lmerTest::lmer(fx, data=subdat)
  #############################################################################
  ## Q3. please help me better understand output of:
  #summary(fm)
  #############################################################################
  #############################################################################
  ## Q4. Goodness of fit. How to asses?
  #qqnorm(resid(fm))
  #qqline(resid(fm))
  #############################################################################
  # compute some model statistics, store in list rho
  rho <- list()
  rho$protein <- protein
  rho$formula <- fx
  rho$model <- fm
  rho$data <- subdat
  # here we compute the unmoderated statistics
  rho$df.prior <- 0
  rho$s2.prior <- 0
  # check if singular (boundary) fit
  rho$isSingular <- lme4::isSingular(fm)
  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$thopt <- lme4::getME(fm, "theta")
  # calculate degrees of freedom and sigma^2:
  av <- anova(fm)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"
  # calcuate symtoptic var-covar matrix
  rho$A <- MSstatsTMT::calcApvarcovar(fx, subdat, rho$thopt,rho$sigma) 
  # we store some proteins statistics in a list
  stats_list <- list()
  # calculate posterior s2
  s2.post <- MSstatsTMT::calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)
  # compute variance-covariance matrix
  vss <- MSstatsTMT::vcovLThetaLM(fx,subdat)
  #############################################################################
  ## Q5. How to define contrast matrix?
  vec = rho$coeff
  vec[] <- 0
  vec[1] <- -1
  vec[2] <- +1
  contrast_matrix <- vec
  #############################################################################
  varcor <- vss(t(contrast_matrix), thpars=c(rho$thopt, rho$sigma))
  vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
  se2 <- as.numeric(contrast_matrix %*% as.matrix(vcov) %*% contrast_matrix)
  # calculate variance
  vcov.post <- varcor$unscaled.varcor * s2.post
  variance <- as.numeric(contrast_matrix %*% as.matrix(vcov.post) %*% contrast_matrix)
  # calculate degrees of freedom
  # given params theta and sigma from lmer and the Acovar
  # NOTE: careful formatting can break things
  g <- MSstatsTMT::mygradient(function(x) vss(t(contrast_matrix), x)$varcor, c(rho$thopt, rho$sigma))
  denom <- t(g) %*% rho$A %*% g
  # compute df.posterior
  df.post <- 2 * (se2)^2 / denom + rho$df.prior # df.post
  ## Q5. FC seems inflated?
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% rho$coeff)[, 1] # coeff
  t <- FC / sqrt(variance) # the Fold change and the variance
  # compute the p-value
  p <- 2 * pt(-abs(t), df = df.post) # t-statistic and df.post
  # compile results
  rho$stats <- data.frame(protein=protein,contrast="Mutant-Control",
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
  		 Tstatistic=t, SE=sqrt(variance), DF=df.post, isSingular=rho$isSingular)
  return(rho)
} #EOF

## ----------------------------------------------------------------------------

n_cores = 23
BiocParallel::register(BiocParallel::SnowParam(n_cores))

## loop to fit all proteins
prots = unique(as.character(msstats_prot$Protein))

results_list <- foreach(protein = prots) %dopar% {
	suppressMessages({
	  try(lmerTestProtein(protein),silent=TRUE)
	})
} # EOL


# collect results
results_df <- bind_rows(sapply(results_list,"[","stats"))

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

# annotate with gene symbols
idx <- match(results_df$protein,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  				 symbol=gene_map$symbol[idx],
  				 .after="protein")

# adjust pvals 
results_df <- tibble::add_column(results_df, 
			 Padjust=p.adjust(results_df$Pvalue,"BH"),
			 .after="Pvalue")

# sort
results_df <- results_df %>% arrange(Pvalue)

# examine top results
results_df %>% head() %>% knitr::kable()
