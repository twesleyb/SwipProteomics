#' js
#'
#' jacaard similarity of two vectors
#'
#' @param x - first vector of items.
#'
#' @param y - second vector of items.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' js(x, y)
js <- function(x, y) {
  s <- length(intersect(x, y)) / length(union(x, y))
  return(s)
}
