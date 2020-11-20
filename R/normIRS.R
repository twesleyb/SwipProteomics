normIRS <- function(tp, controls, robust = FALSE) {
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("data.table", quietly = TRUE)
  # Perform IRS Normalization.
  # Calculate average of QC samples for each experiment (ExpMean).
  tp <- ungroup(tp)
  tp_list <- tp %>%
    filter(Treatment == controls) %>%
    group_by(Accession, Experiment) %>%
    summarize(
      nQC = length(Intensity),
      ExpMean = mean(Intensity)
    ) %>%
    group_split()
  # For every Protein, calculate the global mean of QC samples.
  # Calculate normalization factors to equalize GlobalMean and
  # ExpMean.
  tp_list <- lapply(tp_list, function(x) {
    x$GlobalMean <- mean(x$ExpMean, na.rm = TRUE)
    x$RobustMean <- exp(mean(log(x$ExpMean), na.rm = TRUE))
    x$NormFactor <- x$GlobalMean / x$ExpMean
    x$RobustNormFactor <- x$RobustMean / x$ExpMean
    x$NormQC <- x$NormFactor * x$ExpMean
    x$RobustNormQC <- x$RobustNormFactor * x$ExpMean
    return(x)
  })
  # Collect the data in a df.
  df <- do.call(rbind, tp_list)
  # Use NormFactor's to normalize protein measurements.
  tp_norm <- left_join(tp, df, by = c("Accession", "Experiment"))
  if (robust == TRUE) {
    # Use Robust mean.
    tp_norm$Intensity <- tp_norm$Intensity * tp_norm$RobustNormFactor
  } else {
    # Use Arithmetic mean.
    tp_norm$Intensity <- tp_norm$Intensity * tp_norm$NormFactor
  }
  tp_norm <- tp_norm %>% dplyr::select(colnames(tp))
  tp_norm <- as.data.table(tp_norm)
  return(tp_norm)
}
