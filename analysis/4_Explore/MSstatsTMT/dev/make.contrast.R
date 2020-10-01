#!/usr/bin/env Rscript

.makeContrast <- function(groups) {
###############################################################
## make contrast matrix for pairwise comparisons
###############################################################
#' @keywords internal
  ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast.matrix <- matrix(rep(0, length(groups) * ncomp),
    ncol = length(groups)
  )
  colnames(contrast.matrix) <- groups
  count <- 0
  contrast.matrix.rownames <- NULL
  for (j in seq_len(length(groups) - 1)) {
    for (k in (j + 1):length(groups)) {
      count <- count + 1
      # save row name
      contrast.matrix.rownames <- c(
        contrast.matrix.rownames,
        paste(groups[j], groups[k], sep = "-")
      )
      # set constrast value
      contrast.matrix[count, groups[j]] <- 1
      contrast.matrix[count, groups[k]] <- -1
    }
  }
  rownames(contrast.matrix) <- contrast.matrix.rownames
  return(contrast.matrix)
}
