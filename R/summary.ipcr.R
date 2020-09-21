#' @title Inidividual Parameter Contribution Regression Summary
#' @description This functions returns the estimated coefficients of the individual parameter contribution (IPC) regression equations as well as the corresponding standard errors, t-values, and p-values.
#' @param object an ipcr object.
#' @param parameter a character vector with the target parameters that are to be inspected.
#' @param method a character string indicating the method to be used; for method = "both" standard and iterated IPC regression is used (default), for method = "standard" only standard IPC coefficients are displayed, and for method = "iterated" only iterated coefficients are shown.
#' @return NULL
#' @examples
#' # Show the IPC regression output
#' summary(object = IPC_reg)
#' @export summary.ipcr

summary.ipcr <- function(object, parameter = NULL, method = "iterated") {

  # Check if IPC regression has been run
  if (is.null(object$regression)) {stop("IPC regression has not been run")}

  # Select method
  if (method %in% c("iterated", "both") & is.null(object$iterated_regression)) {
    method <- "standard"
  }

  # Check arguments
  if (!is.null(parameter)) {
    if(!(parameter %in% object$info$parameters)) {
      stop("'parameter' not found")
    }
  }

  if (method != "standard" & method != "iterated" & method != "both") {
    stop("'method' must be either 'standard', 'iterated', or 'both'")
  }

  # Write header
  output <- "IPC regression summary"

  # Add white space
  output[length(output) + 1] <- ""

  # Add information from print.ipcr
  output[length(output) + 1] <- paste("Model under investigation:", object$info$name)
  output[length(output) + 1] <- paste("Model class:", paste(object$info$class, collapse = ", "))
  if (method %in% c("standard", "iterated")) {
    output[length(output) + 1] <- paste("IPC regression method:", object$info$method)
  }
  if (method == "both") {
    output[length(output) + 1] <- paste("IPC regression method:", "standard & iterated")
  }
  output[length(output) + 1] <- paste("Covariates:", paste(object$info$covariates, collapse = ", "))
  if("convergence_criterion" %in% names(object$info)) {
    output[length(output) + 1] <- paste("Convergence criterion:", object$info$convergence_criterion)
  }
  if("iteratations" %in% names(object$info)) {
    output[length(output) + 1] <- paste("Number of iterations:", object$info$iterations)
  }

  # Add whitespace
  output[length(output) + 1] <- ""

  # Select target parameters
  if (is.null(parameter)) {
    parameter <- names(object$regression)
  }

  # Number of target parameters
  q <- length(parameter)
  param_position <- which(names(object$regression) %in% parameter)
  significant <- FALSE

  # Add output for standard IPC regression estimates
  if (method == "standard" | method == "both") {

    output[length(output) + 1] <- "Standard IPC regression parameter estimates:"

    for (i in 1:q) {
      output[length(output) + 1] <- paste("Parameter:", parameter[i])
      lm_summary <- capture.output(summary(object$regression[[param_position[i]]]))
      start_row <- grep(pattern = "Coefficients:", x = lm_summary, fixed = TRUE) + 1
      if (identical(grep(pattern = "---", x = lm_summary, fixed = TRUE), integer(0))) {
        # no significant coefficents
        end_row <- grep(pattern = "Residual standard error", x = lm_summary, fixed = TRUE) - 2
      } else {
        significant <- TRUE
        end_row <- grep(pattern = "---", x = lm_summary, fixed = TRUE) - 1
      }
      number_row <- end_row - start_row
      output[(length(output) + 1):(length(output) + 1 + number_row)] <- lm_summary[start_row:end_row]
    }
    if (isTRUE(significant)) {
      output[length(output) + 1] <- "---"
      output[length(output) + 1] <- "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
    }
  }

  # Create output for iterated IPC regression estimates
  if (method == "iterated" | method == "both") {

    significant <- FALSE

    if(method == "both") {
      output[length(output) + 1] <- ""
    }

    output[length(output) + 1] <- "Iterated IPC regression parameter estimates:"

    for (i in 1:q) {
      output[length(output) + 1] <- paste("Parameter:", parameter[i])
      lm_summary <- capture.output(summary(object$iterated_regression[[param_position[i]]]))
      start_row <- grep(pattern = "Coefficients:", x = lm_summary, fixed = TRUE) + 1
      if (identical(grep(pattern = "---", x = lm_summary, fixed = TRUE), integer(0))) {
        # no significant coefficents
        end_row <- grep(pattern = "Residual standard error", x = lm_summary, fixed = TRUE) - 2
      } else {
        significant <- TRUE
        end_row <- grep(pattern = "---", x = lm_summary, fixed = TRUE) - 1
      }
      number_row <- end_row - start_row
      output[(length(output) + 1):(length(output) + 1 + number_row)] <- lm_summary[start_row:end_row]
    }
    if (isTRUE(significant)) {
      output[length(output) + 1] <- "---"
      output[length(output) + 1] <- "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ' 1"
    }
  }

  writeLines(output)
}
