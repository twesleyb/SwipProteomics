#' sumProt
#'
#' summarize to protein level
#'
#' @export sumProt
#'
#' @importFrom dplyr %>% ungroup group_by summarize

sumProt <- function(tp) {

  # Sum to protein level
  tp$Intensity[is.na(tp$Intensity)] <- 0
  tp <- dplyr::ungroup(tp)
  proteins <- tp %>%
    dplyr::group_by(
      Experiment, Sample, Channel,
      Treatment, Accession
    ) %>%
    dplyr::summarize(
      Peptides = length(Intensity),
      Intensity = sum(Intensity, na.rm = TRUE)
    )
  proteins$Intensity[proteins$Intensity == 0] <- NA

  return(proteins)
}
