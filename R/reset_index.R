#' reset_index
#'
#' reset partition indices
#'
#' @param partition - a named vector describing the partition of a network.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @examples
#' reset_index(partition)
reset_index <- function(partition) {
  x <- partition
  # Start from 0.
  if (min(x) == 0) {
    v <- c(0:length(unique(x)))
  } else {
    # Start from 1.
    v <- c(1:length(unique(x)))
  }
  # Map old indices to new ones.
  namen <- unique(x)
  names(v) <- namen[order(namen)]
  y <- as.numeric(v[as.character(x)])
  names(y) <- names(x)
  return(y)
}
