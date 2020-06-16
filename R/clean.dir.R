clean.dir <- function(dirtydir){
	# Remove any files in a directory.

	all_files <- sapply(dirtydir,function(x) list.files(full.names=TRUE))

	silence({ lapply(all_files,unlink) })
	#return()
}
