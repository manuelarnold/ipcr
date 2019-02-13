#' @title Print ipcr Object Information
#' @description This function prints the some of the main information of an ipcr object.
#' @param object an ipcr object.
#' @return NULL
#' @examples
#' # Print main information
#' print(x = IPC_reg)
#' @export print.ipcr

print.ipcr <- function(x) {
  output <- paste("Model under investigation:", x$Info$name)
  if (x$Info$covariates[1] == "no covariates specified") {
    output[2] <- "No covariates specificied"
  } else if (x$Info$covariates[1] == "missing data in covariates") {
    output[2] <- "Incomplete covariates"
  } else {
    output[2] <- paste("Covariates:", paste(x$Info$covariates, collapse = ", "))
  }
  if("convergence_criterion" %in% names(x$Info)) {
    output[3] <- paste("Convergence criterion:", x$Info$convergence_criterion)
  }
  if("iterate" %in% names(x$Info)) {
    output[4] <- x$Info$iterate
  }
  writeLines(output)
}
