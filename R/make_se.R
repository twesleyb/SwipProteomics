#' Data.frame to SummarizedExperiment object
#' conversion using an experimental design
#'
#' \code{make_se} creates a SummarizedExperiment object
#' based on two data.frames: the protein table and experimental design.
#'
#' @param proteins_unique Data.frame,
#' Protein table with unique names annotated in the 'name' column
#' (output from \code{\link{make_unique}()}).
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @return A SummarizedExperiment object
#' with log2-transformed values.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#' @export
make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)

  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID

  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }

  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                                  colData = expdesign,
                                  rowData = row_data)
  return(se)
}
