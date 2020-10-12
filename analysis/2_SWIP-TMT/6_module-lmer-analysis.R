#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: module-level analysis

# inputs ----------------------------------------------------------------------
# * msstats_prot.rda


# options ---------------------------------------------------------------------
BF_alpha = 0.05


# prepare the environment -----------------------------------------------------
# load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  #library(parallel)
  #library(doParallel)
})

# local imports
devtools::load_all()

## MSstatsTMT internal functions ----------------------------------------------

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
      idx <- names(contr) %in% attr(terms(call$formula), "term.labels")
      call[["contrasts"]] <- contr[idx]
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
  # calc asymptotic variance covariance matrix of variance parameters sigma 
  # and theta
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
  if (min(eigval) < sqrt(.Machine$double.eps)) { # tol ~sqrt(.Machine$double.eps)
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
        dev <- dev + 
		2 * determinant(envff$pp$RX())$modulus - 
		p * log(2 * pi * sigsq)
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


lmerTestModule <- function(msstats_prot,module,contrast_matrix) {
  subdat <- msstats_prot %>% filter(Module == module)
  # the data cannot contain missing values
  if (any(is.na(subdat))) {
    warning("The data cannot contain missing values.")
    return(NULL)
  }
  # fit the model:
  fm <- suppressMessages({try(lmerTest::lmer(fx, data=subdat),silent=TRUE)})
  # FIXME: catch errors/warnings!
  if (inherits(fm,"try-error")) { 
    warning("try-error")
    return(NULL)
  }
  # compute some model statistics, store in list rho
  rho <- list()
  rho$module <- module
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
  # NOTE: utlilizes MSstatsTMT internal function .calcApvar to do calculation
  rho$A <- calcApvarcovar(fx, subdat, rho$thopt,rho$sigma) # 
  # we store the proteins statistics in a list
  stats_list <- list()
  # calculate posterior s2
  s2.post <- calcPosterior(rho$s2, rho$s2_df, rho$s2.prior, rho$df.prior)
  # compuate variance-covariance matrix
  vss <- vcovLThetaLM(fx,subdat)
  varcor <- vss(t(contrast_matrix), c(rho$thopt, rho$sigma))
  vcov <- varcor$unscaled.varcor * rho$s2 # scaled covariance matrix
  se2 <- as.numeric(t(contrast_matrix) %*% as.matrix(vcov) %*% contrast_matrix)
  # calculate variance
  vcov.post <- varcor$unscaled.varcor * s2.post
  variance <- as.numeric(t(contrast_matrix) %*% as.matrix(vcov.post) %*% contrast_matrix)
  # calculate degrees of freedom
  # given params theta and sigma from lmer and the Acovar
  # NOTE: careful formatting can break things
  g <- mygradient(function(x) vss(t(contrast_matrix), x)$varcor, c(rho$thopt, rho$sigma))
  denom <- t(g) %*% rho$A %*% g
  # compute df.posterior
  df.post <- 2 * (se2)^2 / denom + rho$df.prior # df.post
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% rho$coeff)[, 1] # coeff
  t <- FC / sqrt(variance) # the Fold change and the variance
  # compute the p-value
  p <- 2 * pt(-abs(t), df = df.post) # t-statistic and df.post
  # munge comparison
  comparison <- gsub("Genotype","",
		     paste(rev(names(contrast_matrix)),collapse="-"))
  # compile results
  rho$stats <- data.frame(module=module,contrast=comparison,
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
    		 Tstatistic=t, SE=sqrt(variance), DF=df.post, 
		 isSingular=rho$isSingular)
  return(rho)
} #EOF


## load the data --------------------------------------------------------------

# load the data and graph partition
data(swip)
data(gene_map)
data(partition)
data(msstats_prot) # msstats_prot

# annotate data with module membership
membership <- partition[msstats_prot$Protein]
msstats_prot$Module <- membership

# drop modules with fewer than five proteins
to_drop <- names(which(table(partition) < 5))

# drop any unassigned proteins
filt_prot <- msstats_prot %>% 
	filter(!is.na(Module)) %>% 
	filter(Module %notin% to_drop)

# the lmer formula to each module-level subset of the data:
fx <- Abundance ~ 0 + (1|BioFraction) + (1|Protein) + Genotype

#fx <- Abundance ~ 0 + (1|BioFraction) + (1|Protein) + (1|Module) + Genotype

# contrast matrix for Contrl-Mutant comparison:
cm1 <- setNames(c(-1,1), nm=c("GenotypeControl","GenotypeMutant"))

# all modules
modules <- unique(filt_prot$Module)
modules <- modules[modules!=0]

message("\nTotal number of modules: ",length(modules))

# loop to peform analysis for all modules
results_list <- list()
pbar <- txtProgressBar(max=length(modules),style=3)
for (module in modules){
  results_list[[module]] <- lmerTestModule(msstats_prot, module,
					   contrast_matrix=cm1)
  setTxtProgressBar(pbar,match(module,modules))
} #EOL
close(pbar)

# names
names(results_list) <- paste0("M",modules)

# pvalue correction
results_df <- bind_rows(lapply(results_list,function(x) x$stats))
results_df$FDR <- p.adjust(results_df$Pvalue,method="BH")
results_df$PAdjust <- p.adjust(results_df$Pvalue,method="bonferroni")

# annotate modules with module size
m <- as.character(results_df$module)
module_size <- sapply(split(partition,partition),length)[m]
results_df <- tibble::add_column(results_df,size=module_size,.after="module") %>%
	filter(!isSingular)

# drop singular column
results_df$isSingular <- NULL

# module is M#
results_df$module <- paste0("M",results_df$module)

# examine top results:
results_df %>% filter(PAdjust < BF_alpha) %>% arrange(Pvalue) %>% head(5)
message("Wash4c Module: ", paste0("M",partition[swip]))

# examine top results
message("Number of sig modules (PAdjust<0.05): ", 
	sum(results_df$PAdjust < BF_alpha))

# annotate results with module proteins
module_list <- split(names(partition),partition)[-1] # drop M0
fx <- function(x) { # x = module_list[[1]]
	idx <- match(x,gene_map$uniprot)
	return(paste(gene_map$id[idx],collapse="; "))
}
module_prots <- sapply(module_list,fx)
names(module_prots) <- paste0("M",names(module_prots))
results_df$Proteins <- module_prots[results_df$module]

# save module results
fwrite(results_df,file.path(root,"rdata","Module_Results.csv"))

# save sig modules
sig_modules <- results_df$module[results_df$PAdjust < BF_alpha]
myfile <- file.path(root,"data","sig_modules.rda")
save(sig_modules,file=myfile,version=2)
