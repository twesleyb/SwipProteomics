#' save_adjm_as_rda
#'
#' @export
#'

# Define a function that saves adjm as an edge list 
# data.frame as an rda object.
save_adjm_as_rda <- function(adjm,file) {
	diag(adjm) <- 0 # ZERO is smaller than NA
	adjm[lower.tri(adjm)] <- 0
	edges <- reshape2::melt(adjm)
	save(edges,file=file,version=2)
}
