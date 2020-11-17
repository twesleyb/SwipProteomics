#' lmerTestProtein

#' @import lmerTest

#' @importFrom dplyr %>%

#' @export lmerTestProtein 

lmerTestProtein <- function(protein, fx, msstats_prot, contrasts) {

  getIndex <- function(namen, dm = lme4::fixef(fm)) {
    # a helper function to find column index
    idy <- which(grepl(namen, names(dm)))
  }

  # check input args
  stopifnot(inherits(fx, "formula"))
  stopifnot(inherits(protein, "character"))
  stopifnot(inherits(msstats_prot, "data.frame"))

  # input contrasts should be astats_df matrix, numeric vector, or list of such
  if (inherits(contrasts, "matrix")) {
    # rows of matrix are contrasts
    # rownames of matrix are contrast names (comparisons)
    contrast_list <- unlist(apply(contrasts, 1, list), recursive = FALSE)
    stopifnot(all(sapply(contrast_list, sum) == 0))
  } else if (inherits(contrasts, "numeric")) {
    # the sum of positive and neg coeff in contrast should be 0
    stopifnot(sum(contrasts[contrasts < 0], contrasts[contrasts > 0]) == 0)
    contrast_list <- list(contrasts)
  } else if (inherits(contrasts, "list")) {
    contrast_list <- contrasts
    # comparison names will be generated from pos and neg coefficients
  } else {
    stop("problem parsing input 'contrasts'.")
  }

  # subset the data
  idx <- msstats_prot$Protein == protein
  subdat <- msstats_prot[idx, ]
  if (any(is.na(subdat))) {
    # the data cannot contain missing values
    warning("The data cannot contain missing values.")
    return(NULL)
  }

  # fit the model
  fm <- lmerTest::lmer(fx, data = subdat)

  # evaluate statistical comparisons for all contrasts
  stats_list <- list()

  for (i in seq(contrast_list)) {

    # insure contrast matches names(fixef(fm))
    comparison <- contrast_list[[i]]
    contrast <- comparison[names(sort(sapply(names(comparison), getIndex)))]

    # assess contrast
    test_results <- lmerTestContrast(fm, contrast)

    # collect results
    test_results$Protein <- protein

    if (!is.null(attr(contrast_list, "names"))) {
      # use names of contrasts in contrast_list if they exist
      test_results$Contrast <- names(contrast_list)[i]
    }

    stats_list[[i]] <- test_results
  } # EOL for every comparison

  # compile results
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  # sort cols
  stats_df <- stats_df[, c(
    "Protein", "Contrast", "log2FC", "percentControl",
    "Tstatistic", "Pvalue", "SE", "DF", "isSingular"
  )]

  return(stats_df)
} # EOF
