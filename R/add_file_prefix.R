#' add_file_prefix
#' add a prefix to filename
add_file_prefix <- function(file_path, file_name = NULL, output_dir = NULL, 
			width = 2,date=FALSE) {
  # Profide either a full path to a file,
  # or a directory.
	if (is_valid_path(file_path)) {
		output_dir <- dirname(file_path)
	} else {
		output_dir <- getwd()
	}
  indexed_files <- list.files(output_dir, pattern = ("[0-9]{2,4}_"))
  if (length(indexed_files) == 0) {
    last_file <- 0
  } else {
    last_file <- max(as.numeric(sapply(strsplit(indexed_files, "_"), "[", 1)))
  }
  index <- formatC(last_file + 1, width, format = "d", flag = "0")
  if (date) { 
	  prefix <- paste(index, Sys.Date(), sep = "_")
	  output_file <- file.path(output_dir, 
				   paste(prefix, basename(file_path), sep = "_"))
  } else if (!date) {
	  output_file <- file.path(output_dir, 
				   paste(index, basename(file_path), sep = "_"))
  }
  return(output_file)
}

is_valid_path <- function(file_path) {
	# Does file_path look like a path?
	is_path <- grepl("/",file_path)
	if (is_path) {
		# Check that the path exists.
		path_exists <- dir.exists(dirname(file_name))
		if (!path_exists) { 
			# If it doesn't exist, then stop.
			stop(paste("Path to file",
				   file_name,"does not exist."))
		}
	}
	return(is_path)
}
