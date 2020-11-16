#' ggsaveTable
#'
#' save a table with appropriate dimensions
#'
#' @param mytable a grob object
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/13867325/get-width-of-plot-area-in-ggplot2}
#' @keywords grob size height width
#' @export ggsaveTable

ggsaveTable <- function(mytable, myfile) {

  suppressPackageStartupMessages({
    require(grid)
    require(ggplot2)
    require(cowplot)
    require(gridExtra)
  })

  # Function to get absolute size of a grob in inches.
  # Modified from: Hack-R's solution on Stackoverflow, see refs
  # Get height (h) and width (w) of table in inches.
  f <- tempfile()
  png(f)
  h <- convertHeight(sum(mytable$heights), "in", TRUE)
  w <- convertWidth(sum(mytable$widths), "in", TRUE)
  dev.off()
  unlink(f)

  # save
  plot <- plot_grid(mytable)
  ggsave(myfile, plot, width = 1.05 * w, height = 1.1 * h, units = "in")

} #EOF
