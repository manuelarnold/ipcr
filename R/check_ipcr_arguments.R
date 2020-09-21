#' @title Check the specified arguments of ipcr
#' @description Internal function that checks if the arguments specified in the ipcr function are valid.
#' @param fit a fitted model.
#' @param iterated a logical value.
#' @param conv an integer used as stopping criterion for iterated IPC regression.
#' @param linear a logical value.
#' @param max_it the maximum number of iterations for iterated IPC regressions.
#' @export

check_ipcr_arguments <- function(fit, iterated, conv, max_it, linear,
                                 model_data) {

  # Check model class
  if (all(!(class(fit) %in% c("clm", "coxph", "glm", "hurdle", "lavaan", "lm",
                          "mlm", "mlogit", "MxRAMModel", "nls", "polr", "rlm",
                          "survreg", "zeroinfl")))) {
    stop(paste("objects of class", class(fit), "are not supported"))
  }

  # Check iterated argument
  if (!is.logical(iterated)) {
    stop("'iterated' must be TRUE or FALSE")
  }

  # Check conv argument
  if (!conv > 0) {
    stop("'conv' must be > 0")
  }

  # Check maximum number of iterations
  if (!max_it > 0) {
    stop("'max_it' must be > 0")
  }

  # Check linear argument
  if (!is.logical(linear)) {
    stop("'linear' must be TRUE or FALSE")
  }

  # Check if the model data is complete
  if (!all(complete.cases(model_data))) {
    stop("The data used to estimate the structural equation model muss be
         complete. IPCs cannot be computed.")
  }

}

