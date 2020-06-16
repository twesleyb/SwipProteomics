enrichR <- function(gene_list,database,quiet=TRUE){
	# Wrapper around enrichR.
	# Supress output from cat within enrichr with capture.output().
	# FIXME: When not run interactively pbar prodcuces a bunchh of 
	# little pbars...
	# Load the enrichR library.
	suppressPackageStartupMessages({
		library(enrichR)
	})
	# Check input gene_list.
	if (!(is.list(gene_list))) { 
		msg <- "Input gene_list must be a named list of gene symbols."
		stop(msg)
	}
	# Input gene_list should be gene symbols.
	# FIXME: how to test?
	# Check if provided database is available.
	dbs <- listEnrichrDbs()[["libraryName"]]
	db <- tolower(database)
	idx <- match(db,tolower(dbs))
	if (is.na(idx)) {
		# If not an exact match, then try grep.
		idx <- which(grepl(db,tolower(dbs)))
	}
	if (length(idx) > 1) {
		msg <- c("Multiple matching databases, ", 
			 "please select a single database:\n", 
			 paste(dbs[idx],collapse="\n"))
		stop(msg)
	}
	if (length(idx)==0){
		msg <- paste("No matching database found.",
			     "To see all available EnrichR databases,",
			     "try enrichR::listEnrichrDbs().")
		stop(msg)
	}
	db <- dbs[idx]
	# Status report.
	if (!quiet) { 
		n <- length(gene_list)
		message(paste("Analyzing",n,"gene groups for",
			      db,"enrichment..."))
		pbar <- txtProgressBar(max=n,style=3) 
	}
	# Loop to perform enrichment analysis.
	results <- list()
	for (i in 1:length(gene_list)){
		if (!quiet) { setTxtProgressBar(pbar,i) }
		genes <- gene_list[[i]]
		capture.output({ results[[i]] <- enrichr(genes,db) })
	}
	if (!quiet) { close(pbar) }
	# Collect results.
	results <- unlist(results,recursive=FALSE)
	names(results) <- names(gene_list)
	return(results)
}
