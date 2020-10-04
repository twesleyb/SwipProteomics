#!/usr/bin/env Rscript

#' @importFrom lme4 fixef getME
#' @importFrom stats sigma
#' @keywords internal

#' @export

getRho <- function(fm) {
  # utilized by fitLMER to calculate some statistics from the fit model
  rho <- list()
  # calculate coeff, sigma, and theta:
  rho$coeff <- lme4::fixef(fm)
  rho$sigma <- stats::sigma(fm)
  rho$thopt <- lme4::getME(fm, "theta")
  # calculate degrees of freedom and sigma^2:
  av <- anova(fm)
  rho$s2_df <- av$DenDF
  rho$s2 <- av$"Mean Sq" / av$"F value"
  # calcuate symtoptic var-covar matrix
  rho$model <- fm
  rho$A <- .calcApvar(rho)
  return(rho)
} # EOF


# -----------------------------------------------------------------------------


fitLMER <- function(fx, msstats_prot, protein = NULL) {
  # If protein == NULL then all proteins will be fit
  if (!is.null(protein)) {
    check <- is.character(protein) &
      protein %in% msstats_prot$Protein
    stopifnot(check)
    proteins <- protein
  } else {
    # all proteins
    proteins <- unique(as.character(mssatats_prot$Protein))
  }
  # the data cannot contain missing values
  if (any(is.na(msstats_prot %>% filter(Protein %in% proteins)))) {
    # FIXME: report the number or rows
    # FIXME: shouldn't missing values have been imputed by MSstats prev?
    msstats_prot <- msstats_prot[!is.na(msstats_prot$Abundance), ]
    warning("Missing values were removed from input 'data'.")
  }
  # FIXME: need to remove untestable proteins
  # empty list for output of loop for every protein
  fit_models <- list()
  # FIXME: helper and parallelize!
  for (protein in proteins) {
    subdat <- msstats_prot %>% filter(Protein == protein)
    suppressMessages({
      fm <- lmerTest::lmer(fx, data = subdat)
    })
    # compute model stats
    rho <- getRho(fm)
    # add the protein-level data
    rho$data <- subdat
    # check if singular fit
    rho[["isSingular"]] <- lme4::isSingular(fm)
    # return list of rho for every protein
    fit_models[[protein]] <- rho
  } # EOL

  return(fit_models)
} # EOF


# -----------------------------------------------------------------------------

getContrasts <- function(comp, groups) {
  # generate all potential pairwise contrasts
  # matrix with rows = contrasts, columns = levels(Condition)
  # for a given contrast (row) indices of positive coeff = 1
  # for a given contrast (row) indices of negative coeff = -1
  # for a given contrast (row) else = 0
  pw_contrasts <- .makeContrast(groups)
  # pairwise contrasts fx sorted the levels alpha, so the contrasts above,
  # need to be reversed to cover all cases.
  new_rows <- sapply(strsplit(rownames(pw_contrasts), "-"), function(x) {
    paste(x[2], x[1], sep = "-")
  })
  rev_contrasts <- -1 * pw_contrasts
  rownames(rev_contrasts) <- new_rows
  # subset, keep comp
  all_contrasts <- rbind(pw_contrasts, rev_contrasts)
  contrast_matrix <- all_contrasts[comp, ]
  return(contrast_matrix)
}

calcPosterior <- function(s2, s2_df, s2.prior = 0, df.prior = 0) {
  # A function to compute unmoderated posterior s2.post
  s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
  return(s2.post)
}


# -----------------------------------------------------------------------------

testContrasts <- function(fit_list, contrast_matrix, moderated = FALSE) {
  # tests for each protein:
  proteins <- names(fit_list)
  # compute prior df and s2
  if (!moderated) {
    # the unmoderated t-statistic:
    df.prior <- 0
    s2.prior <- 0
  } else if (moderated) {
    idx <- sapply(fit_list, "[", "s2_df") != 0
    eb_input_s2 <- as.numeric(sapply(fit_list, "[", "s2")[idx])
    eb_input_df <- as.numeric(sapply(fit_list, "[", "s2_df")[idx])
    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    if (is.infinite(eb_fit$df.prior)) {
      df.prior <- 0
      s2.prior <- 0
    } else {
      df.prior <- eb_fit$df.prior
      s2.prior <- eb_fit$var.prior
    }
  } # EIE done compution of prior df and s2
  # Loop to perform tests for every protein
  for (protein in proteins) {
    # get protein model and its data from the list of protein fits
    rho <- fit_list[[protein]]
    fm <- rho$model
    subdat <- rho$data
    # annotate if moderated or unmoderated statistics being performed
    rho$isModerated <- moderated
    if (!exists("contrast_list")) {
      # We only need to generate the contrasts once, applied to all proteins
      contrast_list <- lapply(seq(nrow(contrast_matrix)), function(x) {
        y <- contrast_matrix[x, ]
        cm <- .make.contrast.single(fm, y, subdat)
        return(cm)
      })
      names(contrast_list) <- rownames(contrast_matrix)
    }
    # for the given protein loop to perform tests for every comparison
    for (comp in names(contrast_list)) {
      # we store the proteins statistics in a list
      stats_list <- list()
      # the formatted contrast
      cm <- contrast_list[[comp]]
      stats_list$Comparison <- comp
      stats_list$Protein <- protein
      # compute stuff (from MSstatsTMT)
      s2.post <- calcPosterior(s2 = rho$s2, s2_df = rho$s2_df)
      vss <- .vcovLThetaL(fm)
      varcor <- vss(t(cm), c(rho$thopt, rho$sigma))
      vcov <- varcor$unscaled.varcor * rho$s2
      se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
      # calculate variance
      vcov.post <- varcor$unscaled.varcor * s2.post
      variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
      # calculate degrees of freedom
      g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
      denom <- t(g) %*% rho$A %*% g
      df.post <- 2 * (se2)^2 / denom + df.prior
      # compute fold change and the t-statistic
      FC <- (cm %*% rho$coeff)[, 1]
      t <- FC / sqrt(variance)
      # compute the p-value
      p <- 2 * pt(-abs(t), df = df.post)
      # put it all together:
      stats_list$"Log2 Fold Change" <- as.numeric(FC)
      stats_list$"P-value" <- as.numeric(p)
      stats_list$"SE" <- as.numeric(sqrt(variance))
      stats_list$"DF" <- as.numeric(df.post)
      # add to rho
      rho$stats[[comp]] <- stats_list
    } # EOL for each contrast
    # we are still in loop through proteins!
    # coerce stats list to data.frame
    df <- bind_rows(rho$stats)
    rho$stats <- df
    # update fit_list
    fit_list[[protein]] <- rho # updated with stats df
  } # #EOL for each protein

  # fit_list now contains rho$stats a df with protein-wise statistical results
   return(fit_list)

} # EOF


# -----------------------------------------------------------------------------

adjustPvalues <- function(fit_list,adj_method="BH") {
	df <- bind_rows(sapply(fit_list,"[","stats"))
	results_list <- df %>% group_by(Comparison) %>% group_split()
	results_list <- lapply(results_list, function(df) {
		          qval <- p.adjust(x$"P-value", method=adj_method)
	                  df <- tibble::add_column(x, "P-adjust" = qval, 
						   .after="P-value")
			  return(df) })
	# return results in a nice format updated with p.adjust
	return(results_list)
}
