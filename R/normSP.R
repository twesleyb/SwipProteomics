normSP <- function(tp, pool) {

  # Store a copy of the data.
  tp <- ungroup(tp)
  tp_copy <- tp

  # Group pooled together.
  tp$Group <- as.numeric(grepl(paste(pool, collapse = "|"), tp$Sample))
  tp_list <- tp %>%
    group_by(Accession, Group) %>%
    dplyr::summarize(
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      n = length(Intensity), .groups = "drop"
    ) %>%
    as.data.table() %>%
    arrange(Accession, Group) %>%
    group_by(Accession) %>%
    group_split()

  # Loop to calculate normalization factors.
  new_list <- list()
  for (i in 1:length(tp_list)) {
    x <- tp_list[[i]]
    x$NormFactor <- c(x$Mean_Intensity[2] / x$Mean_Intensity[1], 1)
    x$Norm_Mean_Intensity <- x$Mean_Intensity * x$NormFactor
    new_list[[i]] <- x
  }
  tp_list <- new_list

  # Collect in a df.
  df <- do.call(rbind, tp_list) %>%
    dplyr::select(Accession, Group, NormFactor)

  # Merge with input data.
  tp_norm <- left_join(tp, df, by = c("Accession", "Group"))

  # Perform normalization step.
  tp_norm$Intensity <- tp_norm$Intensity * tp_norm$NormFactor
  tp_norm <- tp_norm %>% dplyr::select(colnames(tp_copy))
  tp_norm <- as.data.table(tp_norm)

  # Return the normalized data.
  return(tp_norm)
}
