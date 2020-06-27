normSL <- function(tp, groupBy=NULL){

	suppressPackageStartupMessages({
		library(dplyr)
		library(data.table)
	})

	# Create an expression that will be evaluated to dynamically 
	# group data based on user provided grouping factors. Default's
	# to Sample -- every replicate.
	tp <- ungroup(as.data.frame(tp))
	cmd <- paste0("group_by(tp, ",paste(groupBy,collapse=", "),")")

	# Calculate column sums for grouped samples.
	# FIXME: if .groups is set to 'drop' to avoid Warning msg, then error?
	data_list <- eval(parse(text=cmd)) %>% 
		summarize(Total=sum(Intensity,na.rm=TRUE)) %>%
		group_split()

	# Calculate normalization factors.
	data_SL <- lapply(data_list,function(x) {
				     x$"Mean(Total)" <- mean(x$Total,na.rm=TRUE)
				     x$NormFactor <- x$"Mean(Total)"/x$Total
				     x$Norm <- x$Total*x$NormFactor
			      return(x) }) %>% bind_rows()

	# Collect normalization factors as named vector.
	SL_factors <- data_SL$NormFactor
	names(SL_factors) <- as.character(data_SL[["Sample"]])

	# Normalize measurements by sample factors.
	tp$Intensity <- tp$Intensity * SL_factors[as.character(tp[["Sample"]])]

	return(as.data.table(tp))
}
