sampleConnectivity <- function(tp, log = TRUE) {
  # Calculate Ki (standardized connectivity).
  # The standardized connectivity (Z.K; Oldham et al.,)
  # is a quantity that describes the overall strength of
  # connections between a given node (sample) and all of the other
  # nodes (samples) in a network/adjm.
  # The total connectivity of a Node (sample) is the sum of all of
  # its connections (colSum).
  suppressPackageStartupMessages({
    library(WGCNA)
    library(data.table)
    library(dplyr)
  })
  # Cast tp into matix.
  dm <- tp %>%
    dcast(Accession ~ Sample, value.var = "Intensity") %>%
    as.matrix(rownames = TRUE)
  # Evaluate correlations between samples.
  silence({
    if (log) {
      cormat <- WGCNA::bicor(log2(dm),
        use = "pairwise.complete.obs"
      )
    } else {
      cormat <- WGCNA::bicor(dm,
        use = "pairwise.complete.obs"
      )
    }
  })
  # Calculate adjm.
  adjm <- 0.5 * cormat^2 + 0.5
  # Normalized ki by maximum.
  ki <- colSums(adjm) - 1
  kmax <- max(ki)
  Ki <- ki / kmax
  Kmean <- mean(Ki)
  Kvar <- var(Ki)
  Z.Ki <- (Ki - Kmean) / sqrt(Kvar)
  return(Z.Ki)
}
