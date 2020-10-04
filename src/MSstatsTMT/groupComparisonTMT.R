#' Finding differentially abundant proteins across conditions
#'
#' Tests for significant changes in protein abundance across conditions
#' based on
#' a family of linear mixed-effects models in TMT experiment.  Experimental
#' design of case-control study (patients are not repeatedly measured) is
#' automatically determined based on proper statistical model.
#'
#' @export
#' @param data Name of the output of \code{\link{proteinSummarization}}
#' function. It should have columns named Protein, Mixture, TechRepMixture,
#' Run, Channel, Condition, BioReplicate, Abundance.
#'
#' @param contrast.matrix Comparison between conditions of interests.
#' 1) default is 'pairwise', which
#' compare all possible pairs between two conditions. 2) Otherwise, users can
#' specify the comparisons of interest. Based on the levels of conditions,
#' specify 1 or -1 to the conditions of interests and 0 otherwise.
#' The levels of conditions are sorted alphabetically.
#'
#' @param moderated TRUE will moderate t-statistic; FALSE (default)
#' uses ordinary t-statistic.
#'
#' @param adj.method - adjusted method for multiple comparison.
#' "BH" is default.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels from
#' protein level data.
#'
#' @param remove_empty_channel TRUE(default) removes 'Empty' channels
#' from protein level data.
#'
#' @return returns a data.frame with result of inference
#'
#' @examples
#' data(input.pd)
#' # use protein.summarization() to get protein abundance data
#' quant.pd.msstats <- proteinSummarization(input.pd,
#'   method = "msstats",
#'   global_norm = TRUE,
#'   reference_norm = TRUE
#' )
#'
#' test.pairwise <- groupComparisonTMT(quant.pd.msstats, moderated = TRUE)
#'
#' # Only compare condition 0.125 and 1
#' levels(quant.pd.msstats$Condition)
#'
#' # Compare condition 1 and 0.125
#' comparison <- matrix(c(-1, 0, 0, 1), nrow = 1)
#'
#' # Set the names of each row
#' row.names(comparison) <- "1-0.125"
#'
#' # Set the column names
#' colnames(comparison) <- c("0.125", "0.5", "0.667", "1")
#' test.contrast <- groupComparisonTMT(
#'   data = quant.pd.msstats,
#'   contrast.matrix = comparison,
#'   moderated = TRUE
#' )
groupComparisonTMT <- function(data,
                               contrast.matrix = "pairwise",
                               moderated = FALSE,
                               adj.method = "BH",
                               remove_norm_channel = TRUE,
                               remove_empty_channel = TRUE) {

  ## check input data
  required.info <- c(
    "Protein", "BioReplicate", "Abundance", "Run", "Channel",
    "Condition", "TechRepMixture", "Mixture"
  )

  if (!all(required.info %in% colnames(data))) {
    missedAnnotation <- which(!(required.info %in% colnames(data)))
    missedAnnotation.comb <- paste(required.info[missedAnnotation], collapse = ", ")
    if (length(missedAnnotation) == 1) {
      stop(paste(
        "Please check the required input. ** columns :",
        missedAnnotation.comb, "is missed."
      ))
    } else {
      stop(paste(
        "Please check the required input. ** columns :",
        missedAnnotation.comb, ", are missed."
      ))
    }
  }

  ## remove 'Empty' column : It should not used for further analysis
  if (remove_empty_channel & is.element("Empty", unique(data$Condition))) {
    data <- data[data$Condition != "Empty", ]
    data$Condition <- factor(data$Condition)
  }

  ## remove 'Norm' column : It should not used for further analysis
  if (remove_norm_channel & is.element("Norm", unique(data$Condition))) {
    data <- data[data$Condition != "Norm", ]
    data$Condition <- factor(data$Condition)
  }

  ## remove the rows with NA intensities
  data <- data[!is.na(data$Abundance), ]

  ## Inference
  result <- .proposed.model(data, moderated, contrast.matrix, adj.method)

  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(result)[colnames(result) == "Comparison"] <- "Label"
  colnames(result)[colnames(result) == "adjusted.pvalue"] <- "adj.pvalue"

  return(result)
}
