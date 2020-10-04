#!/usr/bin/env Rscript

load_fun <- function(funcdir) {
	# load .R files in a funcdir
	fun <- list.files(funcdir, pattern="*.R$", full.names=TRUE)
	n <- length(fun)
	if (n==0) { 
		warning("No R files in 'funcdir'.") 
	} else {
		invisible(sapply(fun,source))
	}
} #EOF
