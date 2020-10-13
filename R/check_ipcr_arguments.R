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
  if (!all(stats::complete.cases(model_data))) {
    stop("The data used to estimate the structural equation model muss be
         complete. IPCs cannot be computed.")
  }

}

