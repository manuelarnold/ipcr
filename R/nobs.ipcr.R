#' @title Extract the Number of Observations from an ipcr Object
#' @description Extract the number of observations from an ipcr object.
#' @param object an ipcr object.
#' @param ... arguments passed to methods.
#' @return A single integer.
#' @export

nobs.ipcr <- function(object, ...) {

  if (method == "standard" & is.null(object$regression)) {
    stop("No IPC regression found")
  }

  stats::nobs(object$regression[[1]])

}
