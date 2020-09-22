#' Create an excel workbook.
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
#' @examples
#' write.excel(data, file = "foo.xlsx")
#' @export

write_excel <- function(data, file, rowNames = FALSE, colNames = TRUE, ...) {
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
  wb <- createWorkbook()
  # Loop to add a worksheets:
  for (i in 1:length(list)) {
    df <- as.data.frame(list[[i]])
    addWorksheet(wb, sheetName = names(list[i]))
    writeData(wb, sheet = i, df, rowNames = rowNames, colNames = colNames, ...)
  }
  # Save workbook.
  saveWorkbook(wb, file, overwrite = TRUE)
}
