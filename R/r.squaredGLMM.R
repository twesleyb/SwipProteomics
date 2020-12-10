#!/usr/bin/env Rscript

# NOTE: from MuMin package utils-misc.R utils-models and r.squaredGLMM.R
# Currently this is a mess, because I just pulled everything that
# r.squaredGLMM.merMod required. Thus, much of the code below is not needed
# and could be drastically simplified most likely.

.MuMInEnv <- new.env(parent = baseenv())

testStart <-
  function(...) {
    p <- c(...)
    res <- length(find.package(p, quiet = TRUE)) == length(p)
    if (res) {
      for (a in p) {
        suppressPackageStartupMessages(
          library(a, character.only = TRUE, quietly = TRUE)
        )
      }
      require(.packageName, character.only = TRUE)
      options(na.action = "na.fail")
    }
    res
  }

warnonce <- function(id, ...) {
  if (!isTRUE(get0(flag <- paste0("warned.", as.character(id)[1L]), .MuMInEnv,
    ifnotfound = FALSE
  ))) {
    assign(flag, TRUE, envir = .MuMInEnv)
    cl <- match.call()
    cl$id <- NULL
    cl[[1L]] <- as.name("warning")
    eval.parent(cl)
  }
}

`cry` <-
  function(Call = NA, Message, ..., warn = FALSE, domain = paste0("R-", .packageName)) {
    if (is.character(Call)) {
      Call <- call(Call[1L], sys.call(-1L)[[1L]])
    } else if (is.numeric(Call)) {
      Call <- sys.call(Call - 1L)
    } else if (!is.call(Call) && !is.null(Call)) {
      Call <- sys.call(-1L)
    }
    if (warn) {
      warning(simpleWarning(gettextf(Message, ..., domain = domain), Call))
    } else {
      stop(simpleError(gettextf(Message, ..., domain = domain), Call))
    }
  }


`getElement` <- function(object, name) {
  if (isS4(object)) {
    if (.hasSlot(object, name)) slot(object, name) else NULL
  } else {
    object[[name, exact = TRUE]]
  }
}

# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
  function(x) {
    dfnames <- unlist(lapply(x, colnames), use.names = FALSE)
    uq <- !duplicated(dfnames)
    res <- do.call("cbind", x)[, uq]
    colnames(res) <- dfnames[uq]
    return(res)
  }

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
  function(x) {
    all.colnames <- unique(unlist(lapply(x, colnames), use.names = FALSE))
    x <- lapply(x, function(y) {
      y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
      return(y[all.colnames])
    })
    return(do.call("rbind", x))
  }

`videntical` <-
  function(x) all(vapply(x[-1L], identical, logical(1L), x[[1L]]))

# Check class for each object in a list
`linherits` <- function(x, whats) {
  as.logical(vapply(x, inherits, integer(length(whats)), names(whats),
    which = TRUE
  )) == whats
}

# tries to make a list of element names
`.makeListNames` <- function(x) {
  nm <- names(x)
  lapply(seq_along(x), function(i) {
    if (is.null(nm) || nm[i] == "") {
      switch(mode(x[[i]]),
        call = {
          v <- asChar(x[[i]], width.cutoff = 20L)
          if (length(v) != 1L) v <- sprintf("%s...", v[1L])
          v
        },
        symbol = , name = as.character(x[[i]]),
        NULL = , logical = , numeric = , complex = , character = x[[i]], i
      )
    } else {
      nm[i]
    }
  })
}

# test if dependency chain is satisfied: x[n] can be TRUE only if x[n+1] are also TRUE
`.subset_dc` <- function(...) {
  n <- length(x <- c(...))
  if (n > 1L) all(x[-n] >= x[-1L]) else TRUE
}

# vectorized version of .subset_do (used within subset.model.selection)
`.subset_vdc` <- function(...) apply(cbind(..., deparse.level = 0L), 1L, .subset_dc)

`prettyEnumStr` <- function(x, sep = ", ", sep.last = gettext(" and "), quote = TRUE) {
  n <- length(x)
  if (is.function(quote)) {
    x <- quote(x)
  } else {
    if (identical(quote, TRUE)) quote <- '"'
    if (is.character(quote)) x <- paste0(quote, x, quote)
  }
  paste0(x, if (n > 1L) c(rep(sep, n - 2L), sep.last, "") else NULL,
    collapse = ""
  )
}

# `splitList` <- function (x, k) {
# n <- length(x)
# ret <- unname(split.default(x, findInterval(seq_len(n), seq(0L, n +
# 1L, length = k + 1L))))
# if(k > n) ret <- c(ret, vector(k - n, mode = "list"))
# ret
# }


`.parallelPkgCheck` <- function(quiet = FALSE) {
  # all this is to trick the R-check
  if (!("snow" %in% loadedNamespaces())) {
    if (getRversion() < "2.14.0") {
      if (length(find.package("snow", quiet = TRUE))) {
        do.call("require", list("snow"))
      }
    } else if (length(find.package("parallel", quiet = TRUE))) {
      do.call("require", list("parallel", quiet = TRUE))
    }
  }
  if (!exists("clusterCall", mode = "function")) {
    if (quiet) {
      return(FALSE)
    } else {
      stop("cannot find function 'clusterCall'")
    }
  } else {
    return(TRUE)
  }
}

`clusterVExport` <- local({
  `getv` <- function(obj, env = as.environment(1L)) {
    for (i in names(obj)) assign(i, obj[[i]], envir = env)
  }
  function(cluster, ...) {
    Call <- match.call()
    Call$cluster <- NULL
    Call <- Call[-1L]
    vars <- list(...)
    vnames <- names(vars)
    if (is.null(vnames)) {
      names(vars) <- vapply(Call, asChar, "")
    } else if (any(vnames == "")) {
      names(vars) <- ifelse(vnames == "", vapply(Call, asChar, ""), vnames)
    }
    get("clusterCall")(cluster, getv, vars)
    # clusterCall(cluster, getv, vars)
  }
})

# test if 'x' can be updated (in current environment or on a cluster)
# level is 0/FALSE - no checking, 1 - check if variables and functions exist,
# >1 - reevaluate x and compare with original
`testUpdatedObj` <- function(cluster = NA, x, call = get_call(x),
                             level = 1L, exclude = "subset") {
  if (isTRUE(level)) level <- 2L

  if (level > 0L) {
    xname <- asChar(substitute(x))
    doParallel <- inherits(cluster, "cluster")
    if (doParallel) {
      clusterCall <- get("clusterCall")
      whereStr <- gettext(" in the cluster nodes' environment")
      csapply <- function(...) clusterCall(cluster, "sapply", ...)
    } else {
      whereStr <- ""
      csapply <- function(...) sapply(...)
    }
    if (is.null(call)) stop(gettextf("'%s' has no call component", xname))
    call.orig <- call
    if (!is.null(call$data)) {
      # get rid of formulas, as they are evaluated within 'data'
      call <- call[!sapply(call, function(x) "~" %in% all.names(x))]
      call[exclude] <- NULL
    }

    v <- all.vars(call, functions = FALSE)
    if (!all(z <- unlist(csapply(v, "exists", where = 1L)))) {
      z <- unique(names(z[!z]))
      stop(sprintf(ngettext(
        length(z), "variable %s not found%s",
        "variables %s not found%s"
      ), prettyEnumStr(z, quote = "'"), whereStr))
    }
    vfun <- all.vars(call, functions = TRUE)
    if (!all(z <- unlist(csapply(vfun[!(vfun %in% v)], "exists",
      mode = "function", where = 1L
    )))) {
      zz <- unique(names(z[!z]))
      stop(sprintf(ngettext(
        length(zz), "function %s not found%s",
        "functions %s not found%s"
      ), prettyEnumStr(zz, quote = "'"), whereStr))
    }
    if (level > 1L && !missing(x)) {
      if (doParallel) {
        # XXX: Import: clusterCall
        if (!all(vapply(lapply(clusterCall(cluster, eval, call.orig), all.equal, x), isTRUE, TRUE))) {
          stop(gettextf(
            "'%s' evaluated on the cluster nodes differs from the original one",
            xname
          ))
        }
      } else if (!isTRUE(all.equal(x, update(x)))) {
        stop(gettextf("updated '%s' differ(s) from the original one", xname))
      }
    }
  }
}

`tryCatchWE` <- function(expr) {
  Warnings <- NULL
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      Warnings <<- c(Warnings, list(w))
      invokeRestart("muffleWarning")
    }
  ), warnings = Warnings)
}

# like apply(, 2) but returns a list (does not do any checking)
`applyrns` <- function(X, FUN, ...) {
  n <- nrow(X)
  ret <- vector(n, mode = "list")
  for (i in seq_len(n)) if (!is.null(z <- FUN(X[i, ], ...))) ret[[i]] <- z
  ret
}


## from stats:::format.perc
`format.perc` <-
  function(probs, digits) {
    paste(
      format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
      "%"
    )
  }

return_null <-
  function(...) NULL

## Cheating RCheck:
getFrom <-
  function(pkg, name) {
    get(name, envir = asNamespace(pkg), inherits = FALSE)
  }

# used by 'model.sel' and 'dredge' with argument 'extra'
.get.extras <-
  function(extra, r2nullfit = FALSE) {
    extraExpr <- substitute(extra)
    if (!is.vector(extra)) {
      extraExpr <- call("alist", extraExpr)
      extra <- list(extra)
    }
    if (any(sapply(extra, is.function))) {
      extraExpr[[1L]] <- as.name("alist")
      extra <- eval.parent(extraExpr)
    }
    extraNames <- sapply(extra, function(x) {
      switch(mode(x),
        call = asChar(x[[1L]]), name = asChar(x), character = , x
      )
    })
    if (!is.null(names(extra))) {
      extraNames <- ifelse(names(extra) != "", names(extra), extraNames)
    }
    extra <- structure(as.list(unique(extra)), names = extraNames)
    if (any(i <- vapply(extra, is.language, TRUE))) {
      extra[i] <- lapply(extra[i], eval)
    }

    if (any(c("adjR^2", "R^2") %in% extra)) {
      if (r2nullfit) {
        extra[extra == "R^2"][[1L]] <-
          function(x) {
            r.squaredLR(x,
              null =
                get("nullfit_", parent.frame())
            )
          }
        extra[extra == "adjR^2"][[1L]] <-
          function(x) {
            attr(r.squaredLR(x,
              null =
                get("nullfit_", parent.frame())
            ), "adj.r.squared")
          }
      } else {
        extra[extra == "R^2"][[1L]] <- r.squaredLR
        extra[extra == "adjR^2"][[1L]] <-
          function(x) attr(r.squaredLR(x), "adj.r.squared")
      }
    }
    sapply(extra, match.fun, simplify = FALSE)
  }

## matrix multiplication with option of calculating the diagonal only
matmult <-
  function(x, y, diag.only = FALSE) {
    if (ncol(x) != nrow(y)) stop("non-conformable arguments")
    n1 <- nrow(x)
    n2 <- ncol(y)
    if (diag.only) {
      if (n1 != n2) stop("non-conformable arguments")
      ## >2x faster:
      return(rowSums(x * t(y)))
      # res <- numeric(n1)
      # for(i in seq.int(n1)) res[i] <- sum(x[i, ] * y[, i])
    } else {
      res <- matrix(nrow = n1, ncol = n2)
      for (i in seq.int(n1)) {
        for (j in seq.int(n2)) {
          res[i, j] <- sum(x[i, ] * y[, j])
        }
      }
      return(res)
    }
  }

## matmultdiag(x, ty = y) == diag(x %*% t(y))
matmultdiag <-
  function(x, y, ty = t(y)) {
    if (ncol(x) != ncol(ty)) stop("non-conformable arguments")
    if (nrow(x) != nrow(ty)) stop("result is not a square matrix")
    return(rowSums(x * ty))
  }


tmpvarname <- function(envir, n = 8L) {
  while (exists(x <- paste0(c("*", sample(letters, n), "*"),
    collapse = ""
  ), envir)) {}
  x
}


# Consistent sigma
sigma2 <- function(object) UseMethod("sigma2")
sigma2.default <- function(object) sigma(object)
sigma2.glmmPQL <- function(object) {
  switch(family(object)$family,
    gaussian = , Gamma = object$sigma,
    object$sigma^2
  )
}
sigma2.glmmTMB <- function(object) {
  if (family(object)$family == "nbinom1") sigma(object) + 1 else sigma(object)
}

# VarCorr wrapper returning consistent format (list of named RE variances)
.varcorr <- function(object, ...) UseMethod(".varcorr")
.varcorr.default <- function(object, ...) {
  unclass(lme4::VarCorr(object, ...))
}

# RE model matrix colnames for models with >1 random formulas are prefixed with
# the grouping factor name, e.g. :
# {~ 1 | X1, ~ 1 | X2} has model.matrix columns "X1.(Intercept)", "X2.(Intercept)"
# Need to rename VC matrix dimnames to match them.
.varcorr.lme <- function(object, ...) {
  reStruct <- object$modelStruct$reStruct
  rval <- lapply(reStruct, function(v, sig2) nlme::pdMatrix(v) * sig2, object$sigma^2)
  if ((m <- length(rval)) > 1L) {
    nm <- names(rval)
    for (i in seq.int(m)) {
      dn <- paste(nm[i], dimnames(rval[[i]])[[1L]], sep = ".")
      dimnames(rval[[i]]) <- list(dn, dn)
    }
  }
  rval
}


# Note: currently ONLY FOR CONDITIONAL MODEL
.varcorr.glmmTMB <- function(object, ...) {
  unclass(lme4::VarCorr(object, ...)$cond)
}
.varcorr.glmmadmb <- function(object, ...) {
  suppressWarnings(lme4::VarCorr(object))
}

.numfixef <- function(object, ...) UseMethod(".numfixef")
.numfixef.default <- function(object, ...) lme4::fixef(object, ...)
.numfixef.cpglm <- function(object, ...) coef(object, ...)
.numfixef.glmmTMB <- function(object, ...) lme4::fixef(object, ...)$cond

# RE model matrix
.remodmat <- function(object) UseMethod(".remodmat")

.remodmat.default <-
  function(object) {
    env <- environment(formula(object))
    rval <- lapply(
      .findbars(formula(object)),
      function(f) {
        model.matrix(as.formula(call("~", f[[2L]]), env = env),
          data = model.frame(object)
        )
      }
    )
    rval <- do.call("cbind", rval)
    rval[, !duplicated(colnames(rval)), drop = FALSE]
  }

.remodmat.merMod <-
  function(object) {
    rval <- do.call("cbind", model.matrix(object, type = "randomListRaw"))
    rval[, !duplicated(colnames(rval)), drop = FALSE]
  }

.remodmat.lme <-
  function(object) {
    model.matrix(object$modelStruct$reStruct, data = object$data[rownames(object$fitted), ,
      drop = FALSE
    ])
  }


.nullUpdateWarning <-
  function(message =
             "The null model is correct only if all variables used by the original model remain unchanged.",
           Call = NULL) {
    if (!isTRUE(getOption("MuMIn.noUpdateWarning"))) {
      cry(Call, message, warn = TRUE)
    }
  }

# .nullFitRE: update `object` to intercept only model, keeping original RE terms.
# TODO: reOnlyModelCall or reOnlyFormula
.nullFitRE <- function(object, envir) UseMethod(".nullFitRE")
.nullFitRE.default <-
  function(object, envir = parent.frame()) {
    cl <- getCall(object)
    if (!"formula" %in% names(cl)) {
      stop("cannot create a null model for object without named \"formula\" argument in its call")
    }
    cl$formula <- .nullREForm(formula(object))
    .nullUpdateWarning()
    eval(cl, envir)
  }

.nullFitRE.lme <-
  function(object, envir = parent.frame()) {
    cl <- getCall(object)
    cl$fixed <- update.formula(cl$fixed, . ~ 1)
    if (inherits(object, "glmmPQL")) cl$verbose <- FALSE
    .nullUpdateWarning()
    eval(cl, envir)
  }

# sum up RE variance using VarCorr list
# For RE-intercept identical to a sum of diagonals of VC matrices.
# sum(sapply(lapply(vc, diag), sum))
.varRESum <- function(vc, X) {
  n <- nrow(X)
  sum(sapply(vc, function(sig) {
    mm1 <- X[, rownames(sig), drop = FALSE]
    sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
  }))
}

# update model adding an observation level RE term
.OLREFit <- function(object) UseMethod(".OLREFit")
.OLREFit.default <- function(object) .NotYetImplemented()
.OLREFit.merMod <- function(object) {
  if (!any(sapply(object@flist, nlevels) == nobs(object))) {
    cl <- get_call(object)
    frm <- formula(object)
    nRowData <- eval(call("eval", as.expression(call("NROW", cl$formula[[2L]])),
      envir = cl$data
    ),
    envir = environment(frm),
    enclos = parent.frame()
    )
    fl <- length(frm)
    frx <- . ~ . + 1
    frx[[3L]][[3L]] <- call("(", call("|", 1, call("gl", nRowData, 1)))
    cl$formula <- update.formula(frm, frx)
    object <- tryCatch(eval(cl, envir = environment(frm), enclos = parent.frame()),
      error = function(e) {
        cry(conditionCall(e), conditionMessage(e), warn = TRUE)
        cry(cl, "fitting model with the observation-level random effect term failed. Add the term manually")
      }
    )
    .nullUpdateWarning("The result is correct only if all variables used by the model remain unchanged.")
  }
  object
}


# general function
r2glmm <-
  function(family, vfe, vre, vol, link, pmean, lambda, omega, n) {
    if (inherits(family, "family")) {
      link <- family$link
      family <- family$family
    }

    if (missing(vol) || !is.numeric(vol) || is.na(vol)) {
      vol <- switch(paste(family, link, sep = "."),
        gaussian.identity = vol,
        quasibinomial.logit = ,
        binomial.logit = c(
          theoretical = 3.28986813369645 / n,
          delta = 1 / (n * pmean * (1 - pmean))
        ),
        quasibinomial.probit = ,
        binomial.probit = c(
          theoretical = 1 / n,
          delta =
            6.2831853071795862 / n * pmean * (1 - pmean) *
              exp((qnorm(pmean) / 1.4142135623730951)^2)^2
        ),
        quasibinomial.cloglog = ,
        binomial.cloglog = c(
          theoretical = 1.6449340668482264 / n, #  pi^2 / 6
          delta = pmean / n / log(1 - pmean)^2 / (1 - pmean)
        ),
        Gamma.log = , poisson.log = , quasipoisson.log = , nbinom1.log = c(
          delta = omega / lambda,
          lognormal = log1p(omega / lambda),
          trigamma = trigamma(lambda / omega)
        ),
        quasipoisson.sqrt = , nbinom1.sqrt = ,
        "poisson.mu^0.5" = , poisson.sqrt = c(
          delta = 0.25 * omega
        ),
        nbinom2.log = {
          vdelta <- (1 / lambda) + (1 / omega) # omega is theta
          c(
            delta = vdelta,
            lognormal = log1p(vdelta),
            trigamma = trigamma(1 / vdelta)
          )
        },
        Gamma.inverse = , # c( delta = 1 / nu / lambda^2 ),
        NotImplementedFamily =
          stop("not implemented yet for ", family, " and ", link),
        cry(sys.call(-1L), "do not know how to calculate variance for %s(%s)", family, dQuote(link))
      )
    }

    # print(c(varFE = vfe, varRE = vre, varOL = vol))

    vtot <- sum(vfe, vre)
    matrix(c(vfe, vtot) / (vtot + rep(vol, each = 2L)),
      ncol = 2L, byrow = TRUE, dimnames = list(names(vol), c("R2m", "R2c"))
    )
  }


.binomial.sample.size <-
  function(object) {
    tt <- terms(formula(object))
    y <- model.frame(object)[, rownames(attr(tt, "factors"))[attr(tt, "response")]]
    if (is.null(dim(y))) 1 else mean(rowSums(y))
  }

`r.squaredGLMM` <-
  function(object, null, ...) {
    warnonce(
      "rsquaredGLMM",
      simpleWarning(paste0("'r.squaredGLMM' now calculates a revised statistic. See the help page."))
    )
    UseMethod("r.squaredGLMM")
  }

`r.squaredGLMM.merMod` <-
  function(object, null, envir = parent.frame(), pj2014 = FALSE, ...) {
    if (is.logical(envir)) { # backwards compatibility
      tmp <- envir
      if (!missing(pj2014)) envir <- pj2014
      pj2014 <- tmp
    }

    fam <- family(object)
    # varOL <- lambda <- omega <- NA
    fe <- .numfixef(object)
    ok <- !is.na(fe)
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
    varFE <- var(fitted)

    mmRE <- .remodmat(object)

    ## Note: Argument 'contrasts' can only be specified for fixed effects
    ## contrasts.arg = eval(cl$contrasts, envir = environment(formula(object))))

    vc <- .varcorr(object)

    for (i in seq.int(length(vc))) {
      a <- fixCoefNames(rownames(vc[[i]]))
      dimnames(vc[[i]]) <- list(a, a)
    }
    colnames(mmRE) <- fixCoefNames(colnames(mmRE))

    if (!all(unlist(sapply(vc, rownames), use.names = FALSE) %in% colnames(mmRE))) {
      stop(
        "RE term names do not match those in model matrix. \n",
        "Have 'options(contrasts)' changed since the model was fitted?"
      )
    }

    varRE <- .varRESum(vc, mmRE) # == sum(as.numeric(VarCorr(fm)))

    familyName <- fam$family
    if (substr(familyName, 1L, 17L) == "Negative Binomial") {
      familyName <- "nbinom2"
    }

    if (familyName %in% c(
      "quasipoisson", "poisson", "nbinom1", "nbinom2",
      "binomial", "quasibinomial"
    )) {
      if (missing(null) || !is.object(null)) null <- .nullFitRE(object, envir)
      fixefnull <- unname(.numfixef(null))
    }

    switch(familyName,
      gaussian =
        r2glmm(fam, varFE, varRE, vol = sigma2(object)^2),
      binomial = , quasibinomial = {
        vt <- .varRESum(.varcorr(null), mmRE)
        # XXX: inverse-link seems to give more reasonable value for non-logit
        # links, should inv-logit (plogis) be used here always?
        pmean <- fam$linkinv(fixefnull - 0.5 * vt * tanh(fixefnull * (1 + 2 * exp(-0.5 * vt)) / 6))
        r2glmm(fam, varFE, varRE, pmean = pmean, n = .binomial.sample.size(object))
      }, nbinom2 = {
        lambda <- unname(exp(fixefnull + 0.5 * varRE))
        theta <- if (inherits(object, "glmerMod")) {
          lme4::getME(object, "glmer.nb.theta")
        } else {
          sigma2(object)^-2
        }
        r2glmm(familyName, varFE, varRE, lambda = lambda, omega = theta, link = fam$link)
      }, Gamma = {
        nu <- sigma2(object)^-2
        omega <- 1
        r2glmm(fam, varFE, varRE, lambda = nu, omega = omega)
      }, quasipoisson = , nbinom1 = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        omega <- sigma2(object)^-2
        r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
      }, poisson = {
        vt <- .varRESum(.varcorr(null), mmRE)
        lambda <- unname(exp(fixefnull + 0.5 * vt))
        omega <- 1
        rval <- r2glmm(fam, varFE, varRE, lambda = lambda, omega = omega)
        if (inherits(object, "merMod") &&
          familyName == "poisson" && pj2014) {
          xo <- .OLREFit(object)
          vc <- .varcorr(xo)
          fe <- .numfixef(xo)
          ok <- !is.na(fe)
          fitted <- (model.matrix(xo)[, ok, drop = FALSE] %*% fe[ok])[, 1L]

          n <- nrow(mmRE)
          vname <- names(xo@flist)[sapply(xo@flist, nlevels) == n][1L]
          if (!vname %in% names(vc)) vname <- make.names(vname)
          stopifnot(vname %in% names(vc)) ### !!!
          varresid <- vc[[vname]][1L]
          rval <- rbind(pj2014 = r2glmm(fam, var(fitted),
            .varRESum(vc, mmRE) - varresid,
            vol = log1p(1 / exp(mean(fitted))) + varresid
          )[1L, ], rval)
        }
        rval
      }, r2glmm(fam, varFE, varRE)
    )
  }


`r.squaredGLMM.lme` <-
  function(object, null, ...) r.squaredGLMM.merMod(object, null, ...)

`r.squaredGLMM.glmmTMB` <-
  function(object, null, envir = parent.frame(), ...) {
    fx <- fixef(object) # fixed effect estimates
    if (length(fx$zi) != 0L) { # || length(fx$disp) != 0L)
      stop("r.squaredGLMM cannot (yet) handle 'glmmTMB' object with zero-inflation")
    }
    r.squaredGLMM.merMod(object, null, envir, ...)
  }

`r.squaredGLMM.glmmadmb` <-
  function(object, null, envir = parent.frame(), ...) {
    if (object$zeroInflation) {
      stop("r.squaredGLMM cannot (yet) handle 'glmmADMB' object with zero-inflation")
    }
    r.squaredGLMM.merMod(object, null, envir, ...)
  }

`r.squaredGLMM.lm` <-
  function(object, null, envir = parent.frame(), ...) {
    fam <- family(object)
    ok <- !is.na(coef(object))
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% coef(object)[ok])[, 1L]
    delayedAssign(
      "fixefnull",
      coef(if (missing(null) || !is.object(null)) {
        .nullFitRE(object, envir)
      } else {
        null
      })
    )
    varFE <- var(fitted)
    familyName <- fam$family
    if (substr(familyName, 1L, 17L) == "Negative Binomial") {
      familyName <- "nbinom2"
    }

    switch(familyName,
      gaussian =
        r2glmm(fam, varFE, 0, vol = sigma2(object)^2),
      binomial = , quasibinomial = {
        r2glmm(fam, varFE, 0,
          pmean = fam$linkinv(unname(fixefnull)),
          n = .binomial.sample.size(object)
        )
      }, Gamma = {
        nu <- sigma2(object)^-2
        omega <- 1
        r2glmm(fam, varFE, 0, lambda = nu, omega = omega)
      }, nbinom2 = {
        r2glmm(familyName, varFE, 0, lambda = unname(exp(fixefnull)), omega = sigma2(object)^-2, link = fam$link)
      }, quasipoisson = , nbinom1 = {
        r2glmm(fam, varFE, 0, lambda = unname(exp(fixefnull)), omega = sigma2(object)^2)
      }, poisson = {
        r2glmm(fam, varFE, 0, lambda = unname(exp(fixefnull)), omega = 1)
      }, r2glmm(fam, varFE, 0)
    )
  }


# TODO
`r.squaredGLMM.glmmML` <-
  function(object, null, ...) {
    .NotYetImplemented()
  }

`r.squaredGLMM.cplm` <-
  function(object, null, envir = parent.frame(), ...) {
    fam <- family(object)
    if (!fam$link %in% c("mu^0", "log")) {
      stop("not implemented yet for ", fam$family, " and ", fam$link)
    }

    fe <- .numfixef(object)
    ok <- !is.na(fe)
    fitted <- (model.matrix(object)[, ok, drop = FALSE] %*% fe[ok])[, 1L]
    varFE <- var(fitted)
    if (missing(null) || !is.object(null)) null <- .nullFitRE(object, envir)

    if (inherits(object, "cpglm")) {
      varRE <- vt <- 0
    } else {
      mmRE <- .remodmat(object)
      varRE <- .varRESum(.varcorr(object), mmRE) # == sum(as.numeric(VarCorr(fm)))
      vt <- .varRESum(.varcorr(null), mmRE)
    }
    mu <- unname(exp(.numfixef(null) + 0.5 * vt)) # the same as getting lambda
    phi <- object@phi # the dispersion parameter
    p <- object@p # the index parameter
    varO <- c(delta = phi * mu^(p - 2), lognormal = log1p(phi * mu^(p - 2)))
    r2glmm(NA, varFE, varRE, varO)
  }
fixLogLik <-
  function(ll, object) {
    if (is.null(attr(ll, "nall")) && is.null(attr(ll, "nobs"))) {
      attr(ll, "nobs") <- nobs(object)
    }
    ll
  }

`.getLik` <-
  function(x) {
    if (isGEE(x)) {
      list(logLik = quasiLik, name = "qLik")
    } else {
      list(logLik = logLik, name = "logLik")
    }
  }


`.getRank` <-
  function(rank = NULL, rank.args = NULL, object = NULL, ...) {
    rank.args <- c(rank.args, list(...))

    if (is.null(rank)) {
      x <- NULL # just not to annoy R check
      IC <- as.function(c(alist(x = , do.call("AICc", list(x)))))
      attr(IC, "call") <- call("AICc", as.name("x"))
      class(IC) <- c("function", "rankFunction")
      return(IC)
    } else if (inherits(rank, "rankFunction") && length(rank.args) == 0L) {
      return(rank)
    }

    srank <- substitute(rank, parent.frame())
    if (srank == "rank") srank <- substitute(rank)

    rank <- match.fun(rank)
    ICName <- switch(mode(srank), call = as.name("IC"), character = as.name(srank), name = , srank)
    ICarg <- c(list(as.name("x")), rank.args)
    ICCall <- as.call(c(ICName, ICarg))
    IC <- as.function(c(alist(x = ), list(substitute(
      do.call("rank", ICarg),
      list(ICarg = ICarg)
    ))))

    if (!is.null(object)) {
      test <- IC(object)
      if (!is.numeric(test) || length(test) != 1L) {
        stop("'rank' should return numeric vector of length 1")
      }
    }

    attr(IC, "call") <- ICCall
    class(IC) <- c("function", "rankFunction")
    IC
  }


# sorts alphabetically interaction components in model term names
# if 'peel', tries to remove coefficients wrapped into function-like syntax
# (this is meant mainly for 'unmarkedFit' models with names such as "psi(a:b:c)")
# FIXME: this function is ugly
`fixCoefNames` <-
  function(x, peel = TRUE) {
    if (!length(x)) {
      return(x)
    }
    ox <- x
    ia <- grep(":", x, fixed = TRUE)
    if (!length(ia)) {
      return(structure(x, order = rep.int(1L, length(x))))
    }
    x <- ret <- x[ia]
    if (peel) {
      # case of pscl::hurdle, cf are prefixed with count_|zero_
      if (all(substr(x, 1L, pos <- regexpr("_", x, fixed = TRUE)) %in%
        c("count_", "zero_"))) {
        ret <- substr(ret, pos + 1L, 256L)
        k <- TRUE
        suffix <- ""
      } else { # unmarkedFit with its phi(...), lambda(...) etc...
        k <- grepl("^\\w+\\(.+\\)$", x, perl = TRUE)
        fname <- substring(x[k], 1L, attr(regexpr("^\\w+(?=\\()", x[k],
          perl = TRUE
        ), "match.length"))

        # do not peel off if a function
        k[k] <- !vapply(fname, exists, FALSE, mode = "function", envir = .GlobalEnv)
        if (any(k)) {
          pos <- vapply(x[k], function(z) {
            parens <- lapply(
              lapply(
                c("(", ")"),
                function(s) gregexpr(s, z, fixed = TRUE)[[1L]]
              ),
              function(y) y[y > 0L]
            )
            parseq <- unlist(parens, use.names = FALSE)
            p <- cumsum(rep(c(1L, -1L), sapply(parens, length))[order(parseq)])
            if (any(p[-length(p)] == 0L)) -1L else parseq[1L]
          }, 1L, USE.NAMES = FALSE)
          k[k] <- pos != -1L
          pos <- pos[pos != -1]
          if (any(k)) ret[k] <- substring(x[k], pos + 1L, nchar(x[k]) - 1L)
        }
        suffix <- ")"
      }
    } else {
      k <- FALSE
    }

    ## prepare = replace multiple ':' to avoid splitting by '::' and ':::'
    spl <- expr.split(ret, ":", prepare = function(x) gsub("((?<=:):|:(?=:))", "_", x, perl = TRUE))
    ret <- vapply(lapply(spl, base::sort), paste0, "", collapse = ":")
    if (peel && any(k)) {
      ret[k] <- paste0(substring(x[k], 1L, pos), ret[k], suffix)
    }
    ox[ia] <- ret
    ord <- rep.int(1, length(ox))
    ord[ia] <- sapply(spl, length)
    structure(ox, order = ord)
  }

## like 'strsplit', but ignores split characters within quotes and matched
## parentheses
expr.split <-
  function(x, split = ":",
           paren.open = c("(", "[", "{"), paren.close = c(")", "]", "}"),
           quotes = c("\"", "'", "`"), esc = "\\",
           prepare = NULL) {

    ## error checking:
    # if(length(paren.open) != length(paren.close))
    # 	stop("'paren.open' and 'paren.close' are not of the same length")
    # if(any(test <- vapply(c('paren.open', 'paren.close', 'quotes', 'esc', 'split'), function(x, frame) {
    # 	any(nchar(get(x, frame, inherits = FALSE)) != 1L)
    # }, FALSE, frame = sys.frame()))) {
    # 	stop(sprintf(ngettext(sum(test), "argument %s is not single character",
    # 		"arguments %s are not single character") , prettyEnumStr(names(test)[test])),
    # 		 domain = "R-MuMIn")
    # }

    x0 <- x
    if (is.function(prepare)) x <- prepare(x)
    m <- length(x)
    n <- nchar(x)
    res <- vector("list", m)
    for (k in 1L:m) {
      pos <- integer(0L)
      inquote <- ch <- ""
      inparen <- integer(3L)
      for (i in seq.int(n[k])) {
        chprv <- ch
        ch <- substr(x[k], i, i)
        if (inquote != "") { # in quotes
          if (chprv == esc && ch == esc) {
            ch <- " "
          } else
          if (chprv != esc && ch == inquote) inquote <- ""
        } else {
          inparen[j] <- inparen[j <- (inparen != 0L) & (ch == paren.close)] - 1L
          if (ch %in% quotes) {
            inquote <- ch
          } else if (any(j <- (ch == paren.open))) {
            inparen[j] <- inparen[j] + 1L
          } else if (all(inparen == 0L) && ch == split) {
            pos <- c(pos, i)
          }
        }
      }
      res[[k]] <- substring(x0[k], c(1L, pos + 1L), c(pos - 1L, n[k]))
    }
    res
  }

getResponseFormula <-
  function(f) {
    f <- if (!is.null(tf <- attr(f, "terms"))) {
      formula(tf)
    } else {
      formula(f)
    }
    if ((length(f) == 2L) || (is.call(f[[2L]]) && f[[2L]][[1L]] == "~")) {
      0
    } else {
      f[[2L]]
    }
  }


# Tries to find out whether the models are fitted to the same data
checkIsModelDataIdentical <-
  function(models, error = TRUE) {
    cl <- sys.call(sys.parent())
    err <- if (error) {
      function(x) stop(simpleError(x, cl))
    } else {
      function(x) warning(simpleWarning(x, cl))
    }
    res <- TRUE

    responses <- lapply(models, function(x) getResponseFormula(formula(x)))

    if (!all(vapply(responses[-1L], "==", FALSE, responses[[1L]]))) {
      err("response differs between models")
      res <- FALSE
    }

    # datas <- lapply(models, function(x) get_call(x)$data)
    # XXX: need to compare deparse'd 'datas' due to ..1 bug(?) in which dotted
    #  arguments (..1 etc) passed by lapply are not "identical"
    datas <- vapply(lapply(models, function(x) get_call(x)$data), asChar, "")

    # XXX: when using only 'nobs' - seems to be evaluated first outside of MuMIn
    # namespace which e.g. gives an error in glmmML - the glmmML::nobs method
    # is faulty.
    nresid <- vapply(models, function(x) nobs(x), 1) # , nall=TRUE

    if (!all(sapply(datas[-1L], identical, datas[[1L]])) ||
      !all(nresid == nresid[[1L]])) { # better than 'nresid[-1L] == nresid[[1L]]'
      # XXX: na.action checking here
      err("models are not all fitted to the same data")
      res <- FALSE
    }
    invisible(res)
  }

.checkNaAction <-
  function(x, cl = get_call(x),
           naomi = c("na.omit", "na.exclude"), what = "model") {
    naact <- NA_character_
    msg <- NA_character_

    # handles strings, symbols and calls (let's naively assume no one tries to pass
    # anything else here)
    .getNAActionString <- function(x) {
      if (is.symbol(x)) {
        x <- as.character(x)
      } else if (is.call(x)) {
        x <- eval.parent(x, 2L)
        if (is.symbol(x)) x <- as.character(x)
      }
      return(x)
    }
    # TEST:
    # .checkNaAction(list(call = as.call(alist(fun, na.action = getOption("na.action", default = na.fail)))))
    # .checkNaAction(list(call = as.call(alist(fun, na.action = na.fail))))
    # .checkNaAction(list(call = as.call(alist(fun, na.action = na.omit))))



    if (!is.null(cl$na.action)) {
      naact <- .getNAActionString(cl$na.action)
      if (naact %in% naomi) {
        msg <- sprintf("%s uses 'na.action' = \"%s\"", what, naact)
      }
    } else {
      naact <- formals(eval(cl[[1L]]))$na.action
      if (missing(naact)) {
        naact <- getOption("na.action")
        if (is.function(naact)) {
          statsNs <- getNamespace("stats")
          for (i in naomi) {
            if (identical(get(i, envir = statsNs), naact,
              ignore.environment = TRUE
            )) {
              naact <- i
              break
            }
          }
        }

        naact <- .getNAActionString(naact)
        if (is.character(naact) && (naact %in% naomi)) {
          msg <- sprintf(
            "%s's 'na.action' argument is not set and options('na.action') is \"%s\"",
            what, naact
          )
        }
      } else if (!is.null(naact)) {
        naact <- .getNAActionString(naact)
        if (naact %in% naomi) {
          msg <- sprintf("%s uses the default 'na.action' = \"%s\"", what, naact)
        }
      }
    }
    res <- is.na(msg)
    attr(res, "na.action") <- naact
    attr(res, "message") <- msg
    res
  }

`abbreviateTerms` <-
  function(x, minlength = 4, minwordlen = 1,
           capwords = FALSE, deflate = FALSE) {
    if (!length(x)) {
      return(x)
    }
    if (deflate) {
      dx <-
        # gsub("([\\(,]) *\\w+ *= *(~ *(1 *[\\+\\|]?)?)? *", "\\1", x, perl = TRUE)
        gsub("([,\\(\\[]|^)( *~ *)(1 *([\\|\\+] *)?)?", "\\1",
          gsub("([\\(,]) *\\w+ *= *", "\\1", x, perl = TRUE),
          perl = TRUE
        )
    } else {
      dx <- x
    }

    # .DebugPrint(x)
    s <- strsplit(dx, "(?=[\\W_])", perl = TRUE)
    # remove I(...):
    s <- lapply(s, function(z) {
      z <- if ((n <- length(z)) > 3L && all(z[c(1L, 2L, n)] == c("I", "(", ")"))) {
        z[3L:(n - 1L)]
      } else {
        z
      }
      z[z != " "]
    })
    v <- unique(unlist(s, use.names = FALSE))
    i <- grep("[[:alpha:]]", v, perl = FALSE)
    av <- v
    if (length(i)) {
      tb <- rbindDataFrameList(lapply(s, function(x) {
        as.data.frame(rbind(c(table(x))))
      }))
      tb[is.na(tb)] <- 0L

      if (length(v) > length(i)) {
        minlength <-
          minlength - max(c(0L, apply(
            tb[, v[-i], drop = FALSE], 1L,
            "*", nchar(colnames(tb[, v[-i], drop = FALSE]))
          )))
      }
      n <- min(minlength / rowSums(tb[, v[i], drop = FALSE]))
      if (deflate) {
        repl1 <- c("TRUE" = "T", "FALSE" = "F", "NULL" = "")
        for (j in seq_along(repl1)) av[av == names(repl1)[j]] <- repl1[j]
      }
      av[i] <- abbreviate(av[i], max(n, minwordlen))
      if (capwords) {
        av[i] <- paste0(
          toupper(substring(av[i], 1L, 1L)),
          tolower(substring(av[i], 2L))
        )
      }
    }
    for (j in seq_along(s)) s[[j]] <- paste(av[match(s[[j]], v)], collapse = "")
    names(av) <- v
    structure(unlist(s), names = x, variables = av[i])
  }

`modelDescr` <-
  function(models, withModel = FALSE, withFamily = TRUE,
           withArguments = TRUE, remove.cols = c(
             "formula", "random", "fixed", "model",
             "data", "family", "cluster", "model.parameters"
           )) {
    if (withModel) {
      allTermsList <- lapply(models, function(x) {
        tt <- getAllTerms(x)
        rtt <- attr(tt, "random.terms")
        c(tt, if (!is.null(rtt)) paste0("(", rtt, ")") else NULL)
      })
      allTerms <- unique(unlist(allTermsList))
      abvtt <- abbreviateTerms(allTerms)
      variables <- attr(abvtt, "variables")
      abvtt <- gsub("\\(1 \\| (\\S+)(?: %in%.*)?\\)", "(\\1)", abvtt, perl = TRUE)
      abvtt <- sapply(allTermsList, function(x) {
        paste(abvtt[match(x, allTerms)],
          collapse = "+"
        )
      })
    } else {
      abvtt <- variables <- NULL
    }

    if (withFamily) {
      fam <- sapply(models, function(x) {
        tryCatch(unlist(family(x)[c(
          "family",
          "link"
        )]), error = function(e) character(2L))
      })

      f <- fam[1L, ]
      f[is.na(f)] <- ""
      # f <- vapply(strsplit(f, "(", fixed = TRUE), "[", "", 1L)
      # f[f == "Negative Binomial"] <- "negative.binomial"
      # fam <- cbind(fam, unlist(MASS::negative.binomial(1.345)[c("family", "link")]))
      f <- sub("(?:\\((.*)\\))?$", "(\\1", f)
      f <- paste0(f, ifelse(substring(f, nchar(f)) == "(", "", ","), fam[2, ], ")")
      fam <- f

      # fam[2L, fam[2L, ] ==
      # vapply(unique(f),
      # function(x) {
      # rval <- if(is.na(x)) NA_character_ else formals(get(x))$link[1L]
      # if(!is.character(rval)) NA_character_ else rval
      # }, FUN.VALUE = "")[f]] <- NA_character_
      # j <- !is.na(fam[2L,])
      # fnm <- fam[1L, j]
      # fnm <- ifelse(substring(fnm, nchar(fnm)) != ")",
      # paste0(fnm, "("), paste0(substring(fnm, 1, nchar(fnm) - 1),
      # ", "))
      # fam[1L, j] <- paste0(fnm, fam[2L, j], ")")
    }

    if (withArguments) {
      cl <- lapply(models, get_call)
      haveNoCall <- vapply(cl, is.null, FALSE)
      cl[haveNoCall] <- lapply(cl[haveNoCall], function(x) call("none", formula = NA))
      arg <- lapply(cl, function(x) {
        sapply(x[-1L], function(argval) {
          switch(mode(argval), character = , logical = argval,
            numeric = signif(argval, 3L), asChar(argval)
          )
        })
      })
      arg <- rbindDataFrameList(lapply(lapply(arg, t), as.data.frame))
      arg <- cbind(
        class = as.factor(sapply(lapply(models, class), "[", 1L)),
        arg[, !(colnames(arg) %in% remove.cols), drop = FALSE]
      )
      reml <- rep(NA, length(models))
      if (!is.null(arg$method)) {
        reml <- ((arg$class == "lme" &
          is.na(arg$method)) | arg$method == "REML")
        arg$method <- NULL
      }
      if (!is.null(arg$REML)) reml <- ifelse(is.na(arg$REML), reml, arg$REML == "TRUE")
      arg$REML <- as.factor(reml)

      arg <- as.matrix(arg)
      arg[is.na(arg) | arg == "NULL"] <- ""
      arg <- arg[, apply(arg, 2L, function(x) length(unique(x))) != 1L, drop = FALSE]
      if (ncol(arg)) arg <- gsub("([\"'\\s]+|\\w+ *=)", "", arg, perl = TRUE)
    }
    ret <- as.data.frame(cbind(
      model = abvtt, family = if (withFamily) fam else NULL,
      arg, deparse.level = 0L
    ))
    attr(ret, "variables") <- variables
    ret
  }


family2char <-
  function(x, fam = x$family, link = x$link) {
    if (nchar(fam) > 17L && (substr(fam, 1L, 17) == "Negative Binomial")) {
      theta <- as.numeric(strsplit(fam, "[\\(\\)]")[[1L]][2L])
      paste0("negative.binomial", "(", theta, ",", link, ")")
    } else {
      paste0(fam, "(", link, ")")
    }
  }


`commonCallStr` <-
  function(models, calls = lapply(models, get_call)) {
    x <- lapply(calls, as.list)
    alln <- unique(unlist(lapply(x, names)))
    uniq <- vector("list", length(alln))
    names(uniq) <- alln
    uniq[[1L]] <- lapply(x, "[[", 1L)
    for (i in alln[-1]) uniq[[i]] <- lapply(x, "[[", i)
    uniq <- rapply(uniq, classes = "formula", function(x) {
      environment(x) <- .GlobalEnv
      x
    }, how = "replace")
    uniq <- lapply(uniq, unique)
    nu <- sapply(uniq, length)
    strvarious <- "<*>"

    rval <- lapply(uniq, "[[", 1L)
    j <- sapply(rval, inherits, "formula") & nu > 1L
    for (i in which(j)) {
      response <- getResponseFormula(rval[[i]])
      rval[[i]] <- if (identical(response, 0)) {
        call("~", as.name(sprintf("__%d-rhsform__", nu[i])))
      } else {
        call("~", getResponseFormula(rval[[i]]), as.name(sprintf("__%d-rhsform__", nu[i])))
      }
    }

    notj <- !j & nu > 1
    rval[notj] <- paste0("<", nu[notj], " unique values>")
    if (nu[1L] > 1) rval[[1L]] <- paste(sapply(uniq[[1L]], asChar), collapse = "|")

    rval <- paste(deparse(rval[[1L]], control = NULL),
      "(", paste(names(rval[-1L]), "=", rval[-1L], collapse = ", "), ")",
      sep = ""
    )

    rval <- gsub("`__(\\d+)-rhsform__`", "<\\1 unique rhs>", rval, perl = TRUE)
    rval
  }

updateDeps <-
  function(expr, deps) {
    ret <- list()
    env <- sys.frame(sys.nframe())
    expr <- exprapply0(expr, "dc", function(z) {
      v <- vapply(as.list(z[-1L]), asChar, "")
      n <- length(v)
      k <- match(v, colnames(deps))
      for (i in 2L:n) deps[k[1L:(i - 1L)], k[i]] <- TRUE
      assign("deps", deps, envir = env, inherits = FALSE)
      TRUE
    })
    list(deps = deps, expr = expr)
  }

# tests if smooth terms for variables in gam/gamm models have the same 'k'
testSmoothKConsistency <-
  function(models) {

    # method 2: guess from coefficient names:
    if (inherits(models, "model.selection")) {
      x <- lapply(attr(models, "coefTables"), rownames)
      # XXX: add label 'te(x1,x2)'

      res <- unlist(unname(lapply(x, function(x) {
        s <- grep("^(s|t[ei2])\\(.+\\)\\.\\d+$", x, perl = TRUE)
        if (length(s) != 0L) {
          m <- regexpr("^(?:s|t[ei2])\\((.+)\\)\\.\\d+$", x[s], perl = TRUE)
          cst <- attr(m, "capture.start")[, 1L]
          y <- substring(x[s], cst, cst + attr(m, "capture.length")[, 1L] - 1L)
          tapply(y, y, length)
        } else {
          NULL
        }
      })), recursive = FALSE)
      names(res) <- sapply(lapply(expr.split(names(res), ","), sort),
        paste0,
        collapse = ","
      )
    } else {
      # use information stored in gam objects:
      .getSmoothK <-
        function(x) {
          if (inherits(x, "gamm") || (is.list(x) &&
            (length(x) >= 2L) && identical(names(x)[2L], "gam"))) {
            x <- x[[2L]]
          } else if (!inherits(x, "gam")) {
            return(NULL)
          }
          n <- length(x$smooth)
          rval <- vector("list", n)
          nmv <- character(n)
          for (i in seq_len(n)) {
            y <- x$smooth[[i]]
            if (is.null(y$margin)) {
              nmv[i] <- y$term
              rval[[i]] <- y$bs.dim
            } else {
              nm1 <- vapply(y$margin, `[[`, "", "term")
              o <- order(nm1)
              nmv[i] <- paste0(nm1[o], collapse = ",")
              rval[[i]] <- sapply(y$margin, `[[`, "bs.dim")[o]
            }
            nmv[i] <- paste0(sub("\\(.*", "", y$label), "(", nmv[i], ")")
          }
          names(rval) <- nmv
          rval
        }
      res <- unlist(unname(lapply(models, .getSmoothK)), recursive = FALSE)
    }

    if (!is.null(res)) {
      res <- vapply(split(res, names(res)), function(x) {
        k1 <- x[[1L]]
        for (i in 1L:length(x)) if (!identical(x[[i]], k1)) {
          return(TRUE)
        }
        return(FALSE)
      }, FALSE)


      if (any(res)) {
        warning(
          "smooth term dimensions differ between models for variables ",
          prettyEnumStr(names(res)[res], quote = "'"),
          ". Related coefficients are incomparable."
        )
      }
    }
    invisible()
  }
