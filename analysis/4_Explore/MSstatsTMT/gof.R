#' edgeR::gof

renv::load("~/projects/SwipProteomics")
devtools::load_all(quiet=TRUE)
data(edger_fit)
glmfit = fit

data(swip)
data(msstats_prot)
fx = formula(Abundance ~ 0 + (1|BioFraction) + Genotype)
subdat <- subset(msstats_prot,msstats_prot$Protein==swip)
fm <- lmerTest::lmer(fx, data = subdat)
data(results_list)

## need info from all fits!
stopifnot(is(fit,"lmerModLmerTest") | is(fit, "DGEGLM"))

if (inherits(fit,"lmer")) {
	log_like = stats::logLik(fm)
	deviance = -2*log_like[1]
}


gof <- function (glmfit, pcutoff = 0.1, adjust = "holm", plot = FALSE, 
    main = "qq-plot of residual deviances", ...) {

    require("edgeR", character.only=TRUE, quiet= TRUE)
    require("limma", character.only=TRUE, quiet= TRUE)


    gof.stats <- glmfit$deviance  # numeric length(proteins)
    gof.pvals <- pchisq(gof.stats, df = glmfit$df.residual, lower.tail = FALSE, 
        log.p = FALSE)
    outlier <- p.adjust(gof.pvals, method = adjust) < pcutoff
    if (plot) {
        n <- length(gof.stats)
        z <- zscoreGamma(gof.stats, shape = glmfit$df.residual/2, 
            scale = 2)
        col <- rep_len("black", n)
        col[outlier] <- "blue"
        pch <- rep_len(1, n)
        pch[outlier] <- 16
        qqnorm(z, col = col, pch = pch, main = main, ...)
        abline(0, 1)
    }
    invisible(list(gof.statistics = gof.stats, gof.pvalues = gof.pvals, 
        outlier = outlier, df = glmfit$df.residual[1]))
}

gof <- function (lmerfit, pcutoff = 0.1, plot = FALSE, 
    main = "qq-plot of residual deviances", ...) 
{
    stopifnot(is(glmfit, "lmer"))

    gof.stats <- glmfit$deviance
    #gof.pvals <- pchisq(gof.stats, df = glmfit$df.residual, lower.tail = FALSE, 
    #    log.p = FALSE)
    #outlier <- p.adjust(gof.pvals, method = adjust) < pcutoff
    if (plot) {
        n <- length(gof.stats)
        z <- zscoreGamma(gof.stats, shape = glmfit$df.residual/2, scale = 2)
        col <- rep_len("black", n)
        col[outlier] <- "blue"
        pch <- rep_len(1, n)
        pch[outlier] <- 16
        qqnorm(z, col = col, pch = pch, main = main, ...)
        abline(0, 1)
    }
    invisible(list(gof.statistics = gof.stats, gof.pvalues = gof.pvals, 
        outlier = outlier, df = glmfit$df.residual[1]))
}
