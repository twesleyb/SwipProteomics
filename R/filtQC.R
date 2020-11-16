filtQC <- function(tp, controls = "SPQC", nbins = 5, nSD = 4, quiet = TRUE) {
  # Remove peptides with highly variable QC measurments.
  # Calculate ratios of QC peptides, grouped by Experiment.
  # Imports.
  suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
  })
  ratio_data <- tp %>%
    filter(Treatment == controls) %>%
    group_by(Experiment, Accession, Sequence, Modifications) %>%
    dplyr::summarize(
      Ratio = log2(Intensity[2]) - log2(Intensity[1]),
      Mean = log2(mean(Intensity)),
      N = length(Intensity),
      nMissing = sum(is.na(Intensity)),
      Remove = sum(is.na(Intensity) > 0)
    )
  ratio_data <- ratio_data %>% filter(!Remove)
  # Group ratio data into intensity bins.
  breaks <- quantile(ratio_data$Mean, seq(0, 1, length.out = nbins + 1),
    names = FALSE, na.rm = TRUE
  )
  ratio_data$Bin <- cut(ratio_data$Mean, breaks,
    labels = FALSE, include.lowest = TRUE
  )
  # Summarize intensity bins.
  ratio_df <- ratio_data %>%
    group_by(Bin) %>%
    summarize(
      "Median" = median(Mean),
      "Mean" = mean(Ratio, na.rm = TRUE),
      "Std" = sd(Ratio),
      "N" = sum(!is.na(Ratio)),
      "Min" = mean(Ratio, na.rm = TRUE) - (nSD * sd(Ratio)),
      "Max" = mean(Ratio, na.rm = TRUE) + (nSD * sd(Ratio))
    )
  # Determine if QC measurement is outside percision limits.
  ratio_data$Min <- ratio_df$Min[ratio_data$Bin]
  ratio_data$Max <- ratio_df$Max[ratio_data$Bin]
  out_low <- ratio_data$Ratio < ratio_data$Min
  out_high <- ratio_data$Ratio > ratio_data$Max
  out <- out_low | out_high
  ratio_data$isOutlier <- out
  # Summarize number of outlies per bin.
  nOutliers <- ratio_data %>%
    group_by(Bin) %>%
    summarize(n = sum(isOutlier))
  ratio_df$nOutliers <- nOutliers[["n"]]
  # Collect outlier peptides.
  data_outliers <- ratio_data %>% filter(isOutlier)
  outlier_peptides <- paste(data_outliers$Experiment,
    data_outliers$Accession,
    data_outliers$Sequence,
    data_outliers$Modifications,
    sep = "_"
  )
  # Remove outlier peptides from data.
  ids <- paste(tp$Experiment, tp$Accession,
    tp$Sequence, tp$Modifications,
    sep = "_"
  )
  is_outlier <- ids %in% outlier_peptides
  tp_filt <- tp %>% filter(!is_outlier)
  tp_filt <- as.data.frame(tp_filt)
  # Status report.
  if (!quiet) {
    total_out <- sum(out, na.rm = TRUE)
    message(paste(
      "Total number of", controls,
      "outlier peptides identified:", total_out
    ))
  }
  # Return tidy data.
  return(tp_filt)
}
