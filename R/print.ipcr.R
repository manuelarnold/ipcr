#' @title Print ipcr Object Information
#' @description This function prints the main information of an ipcr object.
#' @param x an ipcr object.
#' @return NULL
#' @export print.ipcr

print.ipcr <- function(x) {
  output <- paste("Model name:", x$info$name)
  output[2] <- paste("Model class:", x$info$class)
  output[3] <- paste("Model parameters:", paste(x$info$parameters, collapse = ", "))
  output[4] <- paste("Covariates:", paste(x$info$covariates, collapse = ", "))
  if(!is.null(x$iterated_regression))
  {
    output[5] <- paste("Convergence criterion:", x$info$conv)
    output[6] <- paste("Iteration:", x$info$iterate)
  }
  writeLines(output)
}
