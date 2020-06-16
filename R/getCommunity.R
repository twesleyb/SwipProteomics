#' getCommunity
#'
#' Get community of proteins that interact with seed proteins.
#'
#' @param graph - igraph object
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
#' getCommunity(graph, seeds)
getCommunity <- function(graph, seeds, k = 2) {
  suppressPackageStartupMessages({
    library(igraph)
    library(dplyr)
  })
  # Check that seeds are in provided graph.
  keep <- seeds %in% names(V(graph))
  if (sum(!keep) > 0) {
    message(paste(
      "Warning:", sum(!keep),
      "of the provided seed genes are not found in the graph.",
      "These will be ignored."
    ))
  }
  seeds <- seeds[keep]
  # Get neighborhood of every seed. Slow if there are many seeds.
  node_list <- adjacent_vertices(graph, seeds)
  adj_nodes <- unlist(sapply(node_list, names))
  # Extract community nodes and create new graph.
  subg <- induced_subgraph(graph, adj_nodes)
  all_nodes <- names(V(subg))
  # Remove nodes with less than 2 connections to seed nodes.
  adjm <- as.matrix(as_adjacency_matrix(subg))
  dm <- adjm[rownames(adjm) %in% seeds, ]
  ds <- apply(dm, 2, sum)
  keep <- names(ds)[ds >= k]
  community <- unique(c(seeds, keep))
  community_graph <- induced_subgraph(graph, community)
  return(community_graph)
}
