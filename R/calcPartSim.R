#' partition_similarity
#'
#' calculate pairwise module assignment similarity statistic
#'
#' @param p1 - partition 1 - named vector
#'
#' @param p2 - partition 2 - named vector
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references Choobdar et al.,2019{https://www.ncbi.nlm.nih.gov/pubmed/31471613}
#'
#' @keywords none
#'
#' @examples
#' partition_similarity(p1, p2)

calcPartSim <- function(p1, p2, ignore = 0, normalized = TRUE) {

  # Calculate similarity of two network partitions, from Choobdar et al., 2019.
  # p1 and p2 are named vectors describing the partitioning of two
  # networks into modules
  # p1 and p2 must contain the same set of nodes

  # Check, partitions should have the same names!
  if (!all(names(p1) %in% names(p2))) {
    stop("Partition vectors should contain identical names.")
  }

  p1[p1 == ignore] <- NA
  p2[p2 == ignore] <- NA

  # Create module assignment df.
  df <- data.table::CJ("NodeA" = names(p1), "NodeB" = names(p1), unique = TRUE)
  df <- df[!df$NodeA == df$NodeB, ] # Remove self-interactions.
  df$NodeA_P1 <- as.numeric(p1[df$NodeA])
  df$NodeB_P1 <- as.numeric(p1[df$NodeB])
  df$NodeA_P2 <- as.numeric(p2[df$NodeA])
  df$NodeB_P2 <- as.numeric(p2[df$NodeB])
  df <- na.omit(df)

  # Check if ProtA and ProtB are in same modules in each partition (Pmk).
  Pmk1 <- as.numeric(p1[df$NodeA] == p1[df$NodeB])
  Pmk2 <- as.numeric(p2[df$NodeA] == p2[df$NodeB])

  # Calculate similarity statistic as:
  # Euclidean inner (dot) product / Product of Euclidean norms.
  dp <- pracma::dot(Pmk1, Pmk2)
  enorm_x <- sqrt(sum(Pmk1^2))
  enorm_y <- sqrt(sum(Pmk2^2))
  ps <- dp / (enorm_x * enorm_y)

  if (normalized) {
    # normalize by percent clustered
    not_clustered <- sum(is.na(p1)) + sum(is.na(p1))
    percent_clust <- 1 - (not_clustered / length(c(p1, p2)))
    ps <- ps * percent_clust
  }
  if (is.na(ps)) {
    stop("Warning: Partition Similarity is NA, returning 0.")
  } else {
    return(ps)
  }
} #EOF
