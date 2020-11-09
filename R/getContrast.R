#' getContrast
#' @export getContrast 

getContrast <- function(fm, pos_coef, neg_coef){

  # create a contrast!
  stopifnot(inherits(fm,"lmerModLmerTest"))
  stopifnot(is.character(pos_coef))
  stopifnot(is.character(neg_coef))

  contrast <- lme4::fixef(fm)

  contrast[] <- 0
  neg_index <- grepl(neg_coef,names(contrast))
  pos_index <- grepl(pos_coef,names(contrast))
  contrast[neg_index] <- -1/sum(neg_index)
  contrast[pos_index] <- +1/sum(pos_index)

  stopifnot(any(contrast > 0))
  stopifnot(any(contrast < 0))
  stopifnot(sum(contrast)==0)

  return(contrast)
}
