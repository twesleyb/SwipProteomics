#' is_valid_path
#' Check if a path to a file exists.
#' @export
is_valid_path <- function(file_path) {
	# Does file_path look like a path?
	is_path <- grepl("/",file_path)
	if (is_path) {
		# Check that the path exists.
		path_exists <- dir.exists(dirname(file_path))
		if (!path_exists) { 
			# If it doesn't exist, then stop.
			stop(paste("Path to file",
				   file_path,"does not exist."))
		}
	}
	return(is_path)
}
