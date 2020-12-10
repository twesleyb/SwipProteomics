#' ggsavePDF
#'
#' function for saving multiple ggplots to single pdf
#'
#' @export ggsavePDF
#'
#' @importFrom ggplotify as.ggplot

ggsavePDF <- function(plots, file) {

  # If not a list, coerce to list.
  if (!inherits(plots, "list")) {
    plot_list <- list(plots)
  } else {
    plot_list <- plots
  }

  # Loop through list, save plots to pdf.
  pdf(file, onefile = TRUE)

  for (i in 1:length(plot_list)) {
    # if ggplot, then use print()
    plot <- plot_list[[i]]
    if (inherits(plot, "ggplot")) {
      print(plot)
      # if else, coerce to ggplot and print
    } else {
      print(ggplotify::as.ggplot(plot))
    }
  }
  invisible(dev.off())
}
