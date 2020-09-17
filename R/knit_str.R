#' knit_str
#'
#' knitr::kable(str(robj)) 
#'
#' @return knits the output of str
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references [StackOverflow](https://stackoverflow.com/questions/44200394)
#'
#' @export
#'
#' @import knitr

knit_str <- function(robj,...){
	require(knitr)
	df <- data.frame(variable = names(robj),
		 class = sapply(robj, typeof),
		 first_values = sapply(robj, function(x) paste0(head(x),  collapse = ", ")),
		 row.names = NULL)
	knitr::kable(df,...)
}
