#' summarize_prot
#' @import data.table dplyr

summarize_prot <- function(tp) {
  # Summarize to protein level.
  # Imports.
  suppressPackageStartupMessages({
    require(dplyr, quietly = TRUE)
    require(data.table, quietly = TRUE)
  })
  # Sum to protein level.
  tp$Intensity[is.na(tp$Intensity)] <- 0
  tp <- ungroup(tp)
  proteins <- tp %>%
    group_by(
      Experiment, Sample, Channel,
      Treatment, Accession
    ) %>%
    summarize(
      Peptides = length(Intensity),
      Intensity = sum(Intensity, na.rm = TRUE)
    )
  proteins$Intensity[proteins$Intensity == 0] <- NA
  return(proteins)
}
