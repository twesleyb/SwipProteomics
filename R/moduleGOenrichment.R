#' moduleGOenrichment
#'
#' A function to perform GO enrichmnet analysis of protein co-expression
#' modules.
#'
#' @param
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' moduleGOenrichment()
moduleGOenrichment <- function(partition, gene_map, GOcollection, 
			       useBackground = NULL, p.adjust = "FDR",
			       partition.ids = "id",
			       alpha=0.05,...) {

  # Function to perform GO enrichment for all modules in a given partition.
  suppressPackageStartupMessages({
    library(org.Mm.eg.db)
    library(anRichment)
  })

  # Default background is all genes provided in partition
  if (is.null(useBackground)) { useBackground = "given" }

  # Create a matrix of module labels to be passed to anRichment.
  modules <- split(partition, partition)
  classLabels <- sapply(names(modules), function(x) partition == x)
  names(modules) <- paste0("M",names(modules))
  colnames(classLabels) <- names(modules)
  logic <- classLabels == TRUE
  for (i in 1:ncol(classLabels)) {
    col_header <- colnames(classLabels)[i]
    classLabels[logic[, i], i] <- col_header
    classLabels[!logic[, i], i] <- "NA"
  }
  classLabels <- classLabels[!duplicated(rownames(classLabels)), ]
  classLabels <- classLabels[,-which(colnames(classLabels) == "M0")]

  # Map protein ids to to entrez.
  idx <- match(rownames(classLabels), gene_map[[partition.ids]])
  entrez <- gene_map$entrez[idx]
  rownames(classLabels) <- entrez

  # Perform GO enrichment.
  GOenrichment <- enrichmentAnalysis(
    classLabels,
    identifiers = entrez,
    refCollection = GOcollection,
    active = NULL,
    inactive = NULL,
    useBackground = useBackground,
    threshold = alpha,
    thresholdType = p.adjust,
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "NA",
    verbose = 0
  )

  # Collect the results.
  GO_results <- list()
  for (r in 1:length(GOenrichment$setResults)) {
    GO_results[[r]] <- GOenrichment$setResults[[r]]$enrichmentTable
  }
  names(GO_results) <- colnames(classLabels)

  # Return results.
  return(GO_results)
}
