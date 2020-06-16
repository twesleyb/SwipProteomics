#' sortInput
sortInput <- function(x,idx){
	# Will base function x do the same thing?
	# Used to sort input for NetRep.
	is.symmetric <- function(x) {
		check <- dim(x)[1] == dim(x)[2]
		if (length(check) == 0) {
			return(FALSE)
		} else {
			return(check)
		}
	}
	# Symmetric matrix.
	if (is.matrix(x) & is.symmetric(x)) {
	       	x <- x[idx,idx] 
	# Non-symmetric matrix.
	} else if (is.matrix(x) & !is.symmetric(x)) {
		# Sort by rows.
		if (which(dim(x)==length(idx))==1) { 
			x <- x[idx,] 
		# Sort by columns.
		} else if (which(dim(x)==length(idx))==2) { 
			x <- x[,idx] 
		} else {
			warning("Length idx must be equal to a dimension of x.")
		}
	# Sort a vector.
	} else if (is.vector(x)) {
	       	x <- x[idx] 
	} else {
		warning("Please provided a matrix or vector as input!")
	}
	# Return the sorted data.
	return(x)
}
