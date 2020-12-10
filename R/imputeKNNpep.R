#' imputeKNNpep
#'
#' @import impute
#' @importFrom dplyr %>%
#' @importFrom data.table as.data.table

imputeKNNpep <- function(tp, groupBy = NULL, samples_to_ignore = NULL,
                         colID = "Abundance", max_percent_missing = 0.5,
                         quiet = TRUE) {

  # Create an expression that will be evaluated to dynamically
  # group data based on user provided grouping factors.
  tp <- ungroup(tp)
  if (!is.null(groupBy)) {
    cmd <- paste0("group_by(tp, ", paste(groupBy, collapse = ", "), ")")
    # Split data into groups as specified by groupBy.
    tp_list <- eval(parse(text = cmd)) %>% group_split()
    namen <- sapply(tp_list, function(x) unique(x[[groupBy]]))
    names(tp_list) <- namen
  } else {
    tp_list <- list(tp)
  }

  # Loop to perform imputation of each group seperately.
  results <- list()
  for (i in 1:length(tp_list)) {
    # Status report.
    if (!quiet) {
      if (!is.null(names(tp_list))) {
        g <- names(tp_list)[i]
      } else {
        g <- i
      }
      message(paste("Imputing missing values in group:", g))
    }
    tp_input <- tp_list[[i]]

    # Cast the tidy data into a matrix.
    df <- as.data.table(tp_input) %>%
      dcast(Experiment + Accession +
        Sequence + Modifications ~ Sample,
      value.var = "Intensity"
      )
    idy <- grep(colID, colnames(df))
    dm <- df %>% dplyr::select(all_of(idy))
    dm <- as.matrix(dm)

    # Determine maximum allowable number of missing values.
    N <- ncol(dm)
    idy <- grep(samples_to_ignore, colnames(dm))
    n <- N - length(idy)
    limit <- n * max_percent_missing

    # Ignore, Don't impute if any QC missing.
    ignore1 <- apply(dm[, idy], 1, function(x) sum(is.na(x)) > 0)

    # Ignore, Don't impute if there are too many missing values.
    ignore2 <- apply(dm[, -idy], 1, function(x) sum(is.na(x)) > limit)
    rows_to_ignore <- c(1:nrow(dm))[ignore1 | ignore2]
    rows_to_impute <- c(1:nrow(dm))[-rows_to_ignore]

    # Check if there are any missing values.
    n_missing <- sum(is.na(dm[rows_to_impute, ]))
    p_missing <- round(100 * (n_missing / length(dm)), 3)

    if (n_missing == 0) {
      # If no missing values, stop.
      stop("No missing values to impute.")
    } else if (!quiet) {
      # Else, status report.
      message(paste0(
        "... Percent missing values: ",
        p_missing, "% (n=",
        n_missing, ")."
      ))
    }

    # Perform KNN imputing, suppress messages with quiet.
    knn_data <- impute.knn(dm[rows_to_impute, ],
      k = 10, rowmax = 0.5, colmax = 0.8
    )
    
    # Put the data back together again.
    dm[rows_to_ignore, ] <- NA
    dm[rows_to_impute, ] <- knn_data$data

    # Collect meta data.
    meta <- df %>%
      dplyr::select(which(colnames(df) %notin% colnames(dm)))

    # Combine meta data with imputed data,
    # and then melt into  a tidy df.
    tpKNN <- melt(cbind(meta, dm), id.vars = colnames(meta))

    # Fix column names.
    idy <- c(ncol(tpKNN) - 1, ncol(tpKNN))
    colnames(tpKNN)[idy] <- c("Sample", "Intensity")

    # Sample should be character.
    tp_input$Sample <- as.character(tp_input$Sample)
    tpKNN$Sample <- as.character(tpKNN$Sample)

    # Combine with input data, and
    # ensure the order of the columns is correct.
    tp_input$Intensity <- NULL
    idy <- colnames(tpKNN)[-which(colnames(tpKNN) == "Intensity")]
    tpImpute <- left_join(tp_input, tpKNN, by = idy)
    tpImpute <- tpImpute %>%
      dplyr::select(c(
        "Experiment", "Sample", "Channel", "Treatment",
        "Accession", "Sequence", "Modifications", "Intensity"
      ))
    # Collect results in a list.
    results[[i]] <- tpImpute
  } # Ends loop.

  # Return tidy df.
  tpImpute <- as.data.table(do.call(rbind, results))
  tpImpute <- ungroup(tpImpute)

  return(tpImpute)
}
