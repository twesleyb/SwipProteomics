#' getContrast
#'
#' create a contrast between a models coefficients
#'
#' @export getContrast
#'
#' @importFrom lme4 fixef

getContrast <- function(fm, pos_coef, neg_coef) {

  # FIXME: doesnt work for intercept scenarios

  # create a contrast!

  stopifnot(is.character(pos_coef))
  stopifnot(is.character(neg_coef))

  if (inherits(fm, "lm")) {
    contrast <- fm$coefficients
  } else if (inherits(fm, "lmerModLmerTest")) {
    contrast <- lme4::fixef(fm)
  } else {
    stop("Input 'fm' must be a 'lm' or 'lmerModLmerTest' object.")
  }

  contrast[] <- 0

  neg_index <- grepl(neg_coef, names(contrast))
  pos_index <- grepl(pos_coef, names(contrast))

  contrast[neg_index] <- -1 / sum(neg_index)
  contrast[pos_index] <- +1 / sum(pos_index)

  stopifnot(any(contrast > 0))
  stopifnot(any(contrast < 0))
  stopifnot(sum(contrast) == 0)

  return(contrast)
}
