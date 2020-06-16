#' gse
#'
#' A function to perform gense set enrichment analysis given a 
#  list of entrez ids.
#'
#' @param gene_list - list of entrez genes. elements of list coorespond to
#' different groups.
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
#' gse()
gse <- function(gene_list,GeneSetcollection, background = "given", ...) {
	suppressPackageStartupMessages({
		library(org.Mm.eg.db)
	        library(anRichment)
	})
	# Create a matrix of module labels to be passed to anRichment.
	genes <- unique(unlist(gene_list))
	keep <- !is.na(genes)
	genes <- genes[keep]
	classLabels <- sapply(names(gene_list), function(x) {
				      genes %in% gene_list[[x]] 
			      })
	rownames(classLabels) <- genes
	# Convert matrix T/F into column names.
	logic <- classLabels == TRUE
	for (i in 1:ncol(classLabels)) {
		col_header <- colnames(classLabels)[i]
		classLabels[logic[, i], i] <- col_header
		classLabels[!logic[, i], i] <- "NA"
	}
	# Perform enrichment analysis.
	temp_results <- enrichmentAnalysis(
					   classLabels,
					   identifiers = genes,
	                                   refCollection = GeneSetcollection,
					   active = NULL,
					   inactive = NULL,
					   useBackground = background,
					   threshold = 0.05,
					   thresholdType = "Bonferroni",
					   getOverlapEntrez = TRUE,
					   getOverlapSymbols = FALSE,
					   ignoreLabels = "NA",
					   verbose = 0
					   )
	# Collect the results.
	results <- list()
	n <- ncol(classLabels)
	if (n==1) {
		results[[r]] <- temp_results$enrichmentTable
	} else {
		for (r in 1:n) {
			results[[r]] <- temp_results$setResults[[r]]$enrichmentTable
		}
	}
	names(results) <- colnames(classLabels)
	# Return results.
	return(results)
} # Ends function.
