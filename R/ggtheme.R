#' ggtheme
#'
#' set theme for ggplots
#'
#' @export ggtheme

ggtheme <- function() {

  suppressPackageStartupMessages({
    require(ggplot2)
  })

  ggtheme <- theme_gray() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(
        color = "black", size = 11,
        face = "bold", hjust = 0.5,
        family = "Arial"
      ),
      axis.title.x = element_text(
        color = "black",
        size = 11, face = "bold",
        family = "Arial"
      ),
      axis.title.y = element_text(
        color = "black",
        size = 11, face = "bold",
        family = "Arial"
      ),
      axis.text.x = element_text(
        color = "black", size = 11,
        angle = 0, hjust = 1.0,
        family = "Arial"
      )
    )
  theme_set(ggtheme)
}
