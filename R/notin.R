#' notin
#'
#' %notin% - the opposite of %in%
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import dplyr
#'
#' @examples
#' x %notin% y
#' @export

`%notin%` <- Negate(`%in%`)
