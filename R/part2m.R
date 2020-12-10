#' part2m

part2m <- function(partition,module_prefix="M") {
	modules <- split(names(partition),partition)
	names(modules) <- paste0(module_prefix,names(modules))
	return(modules)
}
