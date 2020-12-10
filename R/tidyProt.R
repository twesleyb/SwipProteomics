#' tidyProt
#'
#' description: a function to tidy proteomics data.
#'
#' @export tidyProt
#'
#' @importFrom reshape2 melt
#'
#' @importFrom dplyr %>% filter
#'
#' @importFrom data.table as.data.table

tidyProt <- function(raw_data, id.vars, species = NULL,
                     samples = NULL, summary = FALSE) {

  lquote <- function(string, single = TRUE) {
    # Wrap a string in single or double quotes.
    single_quote <- "'"
    double_quote <- "\""
    if (single) {
      # Single quote.
      return(paste0(single_quote, string, single_quote))
    } else {
      # Double quote.
      return(paste0(double_quote, string, double_quote))
    }
  }

  dt <- as.data.table(raw_data) %>%
    reshape2::melt(
      id.vars = id.vars,
      variable.name = "Sample",
      value.name = "Intensity",
      variable.factor = FALSE
    ) # Don't coerce to factor.

  # Remove proteins that do not coorespond to species of interest.
  idx <- grepl(paste0("OS=", species), dt$Description)
  if (!all(idx)) {
    n_out <- length(unique(dt$Accession[!idx]))
    msg <- paste(
      n_out, "proteins are not from", lquote(species),
      "and will be removed."
    )
    warning(msg, call. = FALSE)
    dt <- dt %>% dplyr::filter(grepl(paste0("OS=", species), Description))
  }

  # Insure zeros are NA.
  dt$Intensity[dt$Intensity == 0] <- NA

  # Return tidy df.
  return(as.data.table(dt))
}
