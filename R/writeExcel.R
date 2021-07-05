#' Create an excel workbook
#'
#' @param data - named list of data frames (or equivalent) to be written to
#' an excel document.
#'
#' @param file - output file name.
#'
#' @param rowNames - should rownames of data be kept?
#'
#' @param colNames - should colnames of data be kept?
#'
#' @param '...' - additional arguments passed to openxlsx::writeData()
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords write excel write.excel document xlsx xls
#'
#' @import openxlsx
#'
#' @export writeExcel

writeExcel <- function(data, file, rowNames = FALSE, colNames = TRUE, ...) {
  if (class(data) == "list") {
    list <- data
  } else {
    # Coerce to list.
    list <- list(data)
  }
  # Insure there are names.
  if (is.null(names(list))) {
    names(list) <- paste("Sheet", c(1:length(list)))
  }
  wb <- openxlsx::createWorkbook()
  # Loop to add a worksheets:
  for (i in 1:length(list)) {
    df <- as.data.frame(list[[i]])
    openxlsx::addWorksheet(wb, sheetName = names(list[i]))
    openxlsx::writeData(wb, sheet = i, df, rowNames = rowNames, colNames = colNames, ...)
  }
  # Save workbook.
 openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}
