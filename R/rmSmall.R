#' rmSmall

rmSmall <- function(partition, min_size=5) {
	# set module membership of small modules to 0
	stopifnot(min(partition) == 1)
	too_small <- as.numeric(names(which(table(partition) < 5)))
	idx <- partition %in% too_small
	partition[idx] <- 0
	return(partition)
} #EOF
