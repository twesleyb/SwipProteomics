#' scale01

scale01 <- function(x, ...) { 
	# scale range to 0-1
	z <- (x - min(x, ...)) / (max(x, ...) - min(x, ...))
	return(z)
} 
