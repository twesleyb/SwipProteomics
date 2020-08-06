#' convert_to_adjm
#'
convert_to_adjm <- function(edges) {
	suppressPackageStartupMessages({ library(data.table) })
	adjm <- edges %>% as.data.table() %>% 
		dcast.data.table(Var1 ~ Var2, value.var = "value") %>% 
		as.matrix(rownames="Var1")
	return(adjm)
}
