#!/usr/bin/env Rscript

#' @export mapID

mapID <- function(id_in = "Washc4", identifiers = "symbol", new_ids = "uniprot") {
  # requires: env@gene_map
  # mapping between gene identifiers:
  # e.g. mapID('Rogdi', 'symbol', 'uniprot')
  # maps Rogdi from gene symbol to uniprot ID
  if (length(id_in) == 1) {
    type <- "one"
  }
  if (length(id_in) > 1) {
    type <- "many"
  }
  idx <- switch(type,
    one = grepl(id_in, gene_map[[identifiers]]),
    many = match(id_in, gene_map[[identifiers]])
  )
  return(gene_map[[new_ids]][idx])
} # EOF
