#' write_excel
#'
#' Utility function to write data to excel. Data can be provided as a named list.
#'
#' @param data (list, matrix, dataframe)
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords write write.excel save workbook xlsx xls
#'
#' @examples
#' write_excel()

write_excel <- function(data, file, rowNames = FALSE, colNames = TRUE, ...) {
  suppressPackageStartupMessages({
    require(openxlsx, quietly = TRUE)
  })
  if (inherits(data, "list")) {
    data_list <- data
  } else {
    data_list <- list(data)
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste("Sheet", c(1:length(data_list)))
  }
  wb <- createWorkbook()
  for (i in 1:length(data_list)) {
    df <- data_list[[i]]
    addWorksheet(wb, sheetName = names(data_list[i]))
    writeData(wb, sheet = i, df, rowNames=rowNames, colNames=colNames)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}
