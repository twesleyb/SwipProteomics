#' col2hex
#'
#' convert a color to its hexadecimal code
#'
#' @param color - character
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
col2hex <- function(color, maxValue = 255) {
  ## NOTE: this function no longer works...??? 11/2/20

  z <- grDevices::col2rgb(color)

  hex <- rgb(z[1], z[2], z[3], maxColorValue = maxValue)

  return(hex)
}
