#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all()

data(swip)
data(partition)
data(gene_map)
data(msstats_prot)
data(msstats_results)

myfile <- file.path(root,"rdata","adjm.rda")
load(myfile) # adjm

myfile <- file.path(root,"rdata","ne_adjm.rda")
load(myfile) # ne_adjm

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

#' getMedoid
#'
#' Find representative branch from groups determined by heirarchical clustering
#'
#' @param adjm adjacency matrix for heirarchical clustering
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @examples
#' getMedoid(adjm, h)
getMedoid <- function(adjm, k = NULL, h = NULL, method = "complete") {
  # The medoid of the group is the branch that is closest to
  # all branches in its group.
  hc <- hclust(as.dist(1 - adjm), method)
  hc_part <- cutree(hc, k, h)
  groups <- split(names(hc_part), hc_part)
  medoids <- sapply(seq(groups), function(x) {
	       idy = idx = groups[[x]]
	       col_sums <- apply(adjm[idx,idy], 2, sum)
	       medoid <- names(which(col_sums == min(col_sums)))
  })
  return(medoids)
}


prots <- getMedoid(adjm,k=5)

plots <- lapply(prots, plot_profile)

adjm[prots,prots]
ne_adjm[prots,prots]
partition[prots]

plots[[1]]

plots[[2]]

plots[[3]]

plots[[4]]

plots[[5]]

partition[prots]
