#' ggsaveTable
#'
#' save a table with appropriate dimensions.
#'
#' @param x (grob) a grob object
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/13867325/get-width-of-plot-area-in-ggplot2}
#' @keywords grob size height width
#'
#' @export
#'
#' @examples
#' ggsaveTable(grob)
ggsaveTable <- function(mytable, myfile) {
  suppressPackageStartupMessages({
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(cowplot)
  })
  # Function to get absolute size of a grob in inches.
  # Modified from: Hack-R's solution on Stackoverflow, see refernces.
  # Get height (h) and width (w) of table in inches.
  f <- tempfile()
  png(f)
  h <- convertHeight(sum(mytable$heights), "in", TRUE)
  w <- convertWidth(sum(mytable$widths), "in", TRUE)
  dev.off()
  unlink(f)
  # Save.
  plot <- plot_grid(mytable)
  ggsave(myfile, plot, width = 1.05 * w, height = 1.1 * h, units = "in")
}
