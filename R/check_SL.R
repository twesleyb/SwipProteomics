#' check_SL
#'
# This function checks the column sums of 2 random samples from each
# experiment.
check_SL <- function(tp) {
	df <- tp %>% group_by(Experiment,Sample) %>% 
		dplyr::summarize("Sum(Intensity)"=sum(Intensity,na.rm=T),
				 .groups="drop")
	# Random sample from each group.
	idx <- sapply(split(rownames(df),df$Experiment),sample,2,simplify=F)
	knitr::kable(df[unlist(idx),])
}
