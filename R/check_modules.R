#' check_modules
#'
#' Check modules for evidence of preservation.
#'
#' @param selfPreservation result returned by NetRep.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' check_modules(selfPreservation)
check_modules <- function(x, strength = "strong", stats = c(1:7), alpha = 0.05) {
  # Collect observed values, nulls, and p.values -> p.adj.
  obs <- x$observed[, stats]
  nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[, stats]
  q <- apply(x$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
  q[is.na(q)] <- 1
  # If testing more than one statistic, consider strong or weak preservation.
  fx <- c("strong" = "all", "weak" = "any")[strength]
  if (length(stats) > 1) {
    sig <- apply(q < alpha, 1, eval(fx))
    greater <- apply(obs > nulls, 1, eval(fx))
    less <- apply(obs < nulls, 1, eval(fx))
  } else {
    # If testing a single statistic...
    sig <- q < alpha
    greater <- obs > nulls
    less <- obs < nulls
  }
  # Define preserved, divergent, and ns modules.
  nModules <- length(x$nVarsPresent)
  v <- rep("ns", nModules)
  v[greater & sig] <- "preserved"
  v[less & sig] <- "divergent"
  names(v) <- names(x$nVarsPresent)
  return(v)
} # Ends function.
