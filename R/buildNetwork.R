#' buildNetwork
#'
#' Given a HitPredict data.frame (PPIs) and a vector of entrez IDs, build a
#' protein-protein interaction graph.
#'
#' @param hitpredict
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' buildNetwork()
buildNetwork <- function(hitpredict, mygenes, taxid,quiet=TRUE) {
  # Imports.
  suppressPackageStartupMessages({
    require(getPPIs)
    require(dplyr)
    require(igraph)
  })
  # Use genes passed by user, coerce to integer.
  mygenes <- as.integer(mygenes)
  # Check for NA entries.
  is_NA <- is.na(mygenes)
  if (sum(is_NA) > 0) {
    message(paste(
      "Warning:", sum(is_NA),
      "missing or NA entries found in mygenes.",
      "These will be omitted."
    ))
    mygenes <- mygenes[!is_NA]
  }
  # Get interactions among genes of interest.
  ppis <- hitpredict %>% filter(osEntrezA %in% mygenes & osEntrezB %in% mygenes)
  # keep relevant columns.
  sif <- ppis %>% dplyr::select(
    osEntrezA, osEntrezB, Source_database,
    Interaction_detection_methods, Confidence_score,
    Methods, Interactor_A_Taxonomy, EntrezA, EntrezB,
    Publications
  )
  # Map entrez IDs to gene symbols.
  entrez <- unique(c(sif$osEntrezA, sif$osEntrezB, mygenes))
  # Get organism specific mapping database.
  annotationDBs <- mappingDBs()
  orgDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 1) == taxid])
  names(orgDB) <- sapply(strsplit(names(orgDB), "\\."), "[", 2)
  suppressPackageStartupMessages({
    eval(parse(text = paste0("require(", orgDB[["database"]], ",quietly=TRUE)")))
  })
  osDB <- eval(parse(text = orgDB[["database"]]))
  # Get gene symbols, suppress messages.
  suppressMessages({
    symbols <- AnnotationDbi::mapIds(osDB,
      keys = as.character(entrez),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  })
  # Check that all nodes (entrez) are mapped to gene symbols.
  not_mapped <- entrez[is.na(symbols)]
  if (sum(is.na(symbols)) != 0) {
    message(paste(
      "Warning: Unable to map", length(not_mapped),
      "Entrez ID(s) to gene symbol(s)."
    ))
  }
  # Remove nodes that were not successfully mapped.
  nodes <- symbols
  names(nodes) <- entrez
  nodes <- na.omit(nodes)
  # Node attribute table.
  noa <- data.table::data.table(entrez = names(nodes), symbol = nodes)
  rownames(noa) <- noa$entrez
  # Insure that any interactions among unmapped genes are removed.
  out <- sif$osEntrezA %in% not_mapped | sif$osEntrezB %in% not_mapped
  sif <- subset(sif, !out)
  # Build igraph object.
  g <- graph_from_data_frame(sif, directed = FALSE, vertices = noa)
  g <- simplify(g) # remove any redundant edges.
  # Change node names to symbols.
  # g <- set.vertex.attribute(g, "name", value = vertex_attr(g, "symbol"))
  # Status report.
  nNodes <- format(length(V(g)), 1, nsmall = 1, big.mark = ",")
  nEdges <- format(length(E(g)), 1, nsmall = 1, big.mark = ",")
  if (! quiet) {
	  message(paste(nEdges, "edges identified among", nNodes, "nodes!"))
  }
  # Return igraph object.
  return(g)
}
