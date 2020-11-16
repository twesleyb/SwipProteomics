#' col2hex
#'
#' Convert a color to its hexadecimal code.
#'
#' @param color - character
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export col2hex
#'
#' @examples
#' col2hex(color, maxValue = 255)
col2hex <- function(color, maxValue = 255) {
  ## NOTE: this function no longer works...??? 11/2/20
  z <- grDevices::col2rgb(color)
  hex <- rgb(z[1], z[2], z[3], maxColorValue = maxValue)
  return(hex)
}
