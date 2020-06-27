imputeQRILC <- function(tidy_prot,ignore="QC",rowmax=0.5,colmax=0.8,
			  quiet=TRUE){

	# Determine how many missing values each protein has,
	# and then determine the permisible number of missing values 
	# for any given protein.

	# Imports.
	suppressPackageStartupMessages({
		library(dplyr)
		library(tibble)
		library(data.table)
		library(imputeLCMD) # For QRILC method.
	})

	# Store a copy of the input data.
	tp <- tp_in <- tidy_prot

	# How many missing values are there?
	tp$Intensity[tp$Intensity==0] <- NA
	N_missing <- sum(is.na(tp$Intensity))
	if (N_missing ==0) { 
		message(paste("There are no missing values.",
			   "Returning untransformed data."))
		return(tp_in)
	}

	# Separate data to be imputed.
	if (!is.null(ignore)) {
		# Do not impute:
		tp_ignore <- tp %>% filter(grepl(ignore,Sample)) %>% 
			dplyr::select(Accession,Sample,Intensity) %>% 
			as.data.table
		# Data to be imputed:
		tp_impute <- tp %>% filter(!grepl(ignore,Sample)) %>% 
			as.data.table
	} else {
		tp_ignore <- NULL
	}

	# Cast the data into a matrix.
	dm <- tp_impute %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)

	# Don't impute rows (proteins) with too many missing values.
	n_missing <- apply(dm,1,function(x) sum(is.na(x)))
	limit <- ncol(dm) * rowmax
	rows_to_ignore <- n_missing > limit
	n_ignore <- sum(rows_to_ignore)
	if (n_ignore > 0) {
		msg <- paste(n_ignore,"proteins have more than",limit,
			     "missing values, these values cannot be imputed,",
			     "and will be ignored.")
		warning(msg,call.=FALSE)
	}
	# Total number of missing values to be imputed.
	n_imputed <- sum(is.na(dm[!rows_to_ignore,]))
	if (!quiet){
		message(paste("There are",n_imputed, "missing",
		      "values that will be replaced by imputing."))
	}

	# Perform imputing.
		# 
	if (quiet) {
		# Suppress output from impute.knn.
		silence({
			data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
					       k,rowmax,colmax)
		})

# FIXME: tests
data_knn <- impute.QRILC(log2(dm[!rows_to_ignore,]))

	} else {
		data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
				       k,rowmax,colmax)
	}

	# Collect the imputed data.
	dm_knn <- dm
	dm_knn[!rows_to_ignore,] <- 2^data_knn$data

	# Melt into tidy df.
	dt_knn <- as.data.table(dm_knn,keep.rownames="Accession")
	tp_imputed <- melt(dt_knn,id.vars="Accession",
			  variable.name="Sample",value.name="Intensity")

	# Combine with any samples that were ignored.
	tp_imputed <- rbind(tp_ignore,tp_imputed)

	# Combine with input meta data.
	tp_in$Intensity <- NULL
	tp_in$Sample <- as.character(tp_in$Sample)
	tp_imputed$Sample <- as.character(tp_imputed$Sample)
	tp_out <- left_join(tp_in,tp_imputed,by=c("Sample","Accession")) %>%
		as.data.table

	return(tp_out)
}
