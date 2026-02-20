check_arguments_ipcr <- function(fit, predictors, analytic) {

  # Check argument fit ----

  ## Define accepted model classes
  accepted_classes <- c("clm", "coxph", "glm", "hurdle", "lavaan", "lm",
                        "merMod", "mlm", "mlogit", "MxModel", "MxRAMModel",
                        "nls", "polr", "rlm", "survreg", "zeroinfl")

  ## Check if the object has one of the accepted classes
  if (!inherits(fit, what = accepted_classes)) {
    stop(paste("objects of class", class(fit), "are not supported."))
  }

  if (!is.vector(predictors) && !is.matrix(predictors) &&
       !is.data.frame(predictors)) {
    stop("Argument 'predictors' must be a vector, matrix, or data.frame.")
  }


  # Check argument predictors ----

  ## Ensure predictors have at least one column (if matrix or data.frame)
  if (is.matrix(predictors) || is.data.frame(predictors)) {
    if (ncol(predictors) == 0) {
      stop("Argument 'predictors' must contain at least one column.")
    }
  }

  ## Check for missing values
  if (any(is.na(predictors))) {
    warning("Argument 'predictors' contains missing values. Incomplete rows of
            the data will not be used. Consider handling the missing values
            before running ipcr.")
  }

  ## Ensure categorical variables are properly encoded (only relevant for data
  ## frames)
  if (is.data.frame(predictors)) {
    non_numeric_cols <- sapply(predictors, function(col) !is.numeric(col) &&
                                 !is.factor(col))
    if (any(non_numeric_cols)) {
      stop("Argument 'predictors' contains character variables. Convert them to
           factors or dummy variables.")
    }
  }


  # Check argument analytic ----

  ## Check if analytic argument
  if (!is.logical(analytic)) {
    stop("'analytic' must be TRUE or FALSE")
  }

}

check_arguments_ipcr_it <- function(fit, predictors, analytic, conv,
                                    learning_rate, max_it, iteration_info) {
  if (is.null(predictors)) {
    stop("No predictors were provided.")
  }
}


check_arguments_ipcr_it_MxModel <- function(fit, predictors, analytic, conv,
                                            learning_rate, max_it,
                                            iteration_info) {

  if (!inherits(fit, "MxModel"))
    stop("`fit` must be an OpenMx MxModel (or subclass).")

  if (is.null(fit$data$observed) | fit$data$type != "raw") {
    stop("Model must contain raw data (mxData(..., type='raw')).")
  }

}

check_arguments_ipcr_it_MxRAMModel <- function() {

}
