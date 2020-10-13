#' @export

print.ipcr <- function(x, ...) {
  writeLines(x$output$info)
}
