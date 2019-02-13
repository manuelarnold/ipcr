#' @title Check the specified arguments of ipcr
#' @description Internal function that checks if the arguments specified in the ipcr function are valid.
#' @param fit a structural equation model estimated with OpenMx.
#' @param covariates a data frame or matrix with covariates.
#' @param iterated a logical value.
#' @param conv an integer used as stopping criterion for iterated IPC regression.
#' @param max_it the maximum number of iterations for iterated IPC regressions.
#' @export

get_arguments_checked <- function(fit, covariates, iterated, conv, max_it) {

  if (class(fit) != "MxRAMModel") {
    stop("'fit' must be RAM-type OpenMx model object")
  }

  if (!is.null(covariates) & try(class(as.matrix(covariates)), silent = T) != "matrix") {
    stop("'covariates' must be a data frame or matrix")
  }

  if (!is.logical(iterated)) {
    stop("'iterated' must be TRUE or FALSE")
  }

  if (!conv > 0) {
    stop("'conv' must be > 0")
  }

  if (!max_it > 0) {
    stop("'max_it' must be > 0")
  }

  if (!all(complete.cases(fit$data$observed[, fit$manifestVars, drop = FALSE]))) {
    stop("The data used to estimate the structural equation model muss be
         complete. IPCs cannot be computed.")
  }

}

