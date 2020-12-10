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

write_excel <- function(mydata, myfile, rowNames = FALSE, colNames = TRUE, ...) {
  # imports
  suppressPackageStartupMessages({
    require(openxlsx)
  })
  # check that input data is a list
  if ("list" %in% class(mydata)) {
    mylist <- mydata
  } else {
    # Coerce to list
    mylist <- list(mydata)
  }
  # Insure there are names.
  if (is.null(names(mylist))) {
    names(mylist) <- paste("Sheet", c(1:length(mylist)))
  }
  wb <- openxlsx::createWorkbook()
  # Loop to add a worksheets:
  for (i in 1:length(mylist)) {
    df <- as.data.frame(mylist[[i]])
    openxlsx::addWorksheet(wb, sheetName = names(mylist[i]))
    openxlsx::writeData(wb,
      sheet = i, df,
      rowNames = rowNames, colNames = colNames, ...
    )
  }
  # Save workbook.
  openxlsx::saveWorkbook(wb, myfile, overwrite = TRUE)
}
