#' @title Inidividual Parameter Contribution Regression Summary
#' @description This functions returns the estimated coefficients of the individual parameter contribution (IPC) regression equations as well as the corresponding standard errors, t-values, and p-values.
#' @param object an ipcr object.
#' @param parameter a character vector with the target parameters that are to be inspected.
#' @param method a character string indicating the method to be used; for method = "both" initial and iterated IPC regression is used (default), for method = "initial" only initial IPC coefficients are displayed, and for method = "iterated" only iterated coefficients are shown.
#' @return NULL
#' @examples
#' # Show the IPC regression output
#' summary(object = IPC_reg)
#' @export summary.ipcr

summary.ipcr <- function(object, parameter = NULL, method = "both") {

  # Check arguments
  if (!is.null(parameter)) {
    if(!(parameter %in% object$Info$parameters)) {
      stop("'parameter' not found")
    }
  }

  if (method != "initial" & method != "iterated" & method != "both") {
    stop("'method' must be either 'initial', 'iterated', or 'both'")
  }

  # Write header
  output <- "Individual parameter contribution regression summary"

  # Add white space
  output[length(output) + 1] <- ""

  # Add information from print.ipcr
  output[length(output) + 1] <- paste("Model under investigation:", object$Info$name)
  if (object$Info$covariates[1] == "no covariates specified") {
    output[length(output) + 1] <- "No covariates specificied"
  } else if (object$Info$covariates[1] == "missing data in covariates") {
    output[length(output) + 1] <- "Incomplete covariates"
  } else {
    output[length(output) + 1] <- paste("Covariates:", paste(object$Info$covariates, collapse = ", "))
  }
  if("convergence_criterion" %in% names(object$Info)) {
    output[length(output) + 1] <- paste("Convergence criterion:", object$Info$convergence_criterion)
  }
  if("iterate" %in% names(object$Info)) {
    output[length(output) + 1] <- object$Info$iterate
  }

  # Add whitespace
  output[length(output) + 1] <- ""

  # Check if initial IPC regression parameters have been estimated
  if (method == "initial") {
    if ("InitialIPCRegression" %in% names(object)) {
      method_checked <- "initial"
    } else {
      output[length(output) + 1] <- warning("Error: Initial IPC regression paramters have not been estimated.")
    }
  }

  # Check if iterated IPC regression parameters have been estimated
  if (method == "iterated") {
    if ("IteratedIPCRegression" %in% names(object)) {
      method_checked <- "iterated"
    } else {
      output[length(output) + 1] <- warning("Error: Iterated IPC regression paramters have not been estimated.")
    }
  }

  # Check if both initial and iterated IPC regression parameters have been estimated
  if (method == "both") {
    if (all(c("InitialIPCRegression", "IteratedIPCRegression") %in% names(object))) {
      method_checked <- "both"
    } else {
      if (!("InitialIPCRegression" %in% names(object))) {
        output[length(output) + 1] <- warning("Error: Initial IPC regression paramters have not been estimated.")
      }
      if (!("IteratedIPCRegression" %in% names(object))) {
        output[length(output) + 1] <- warning("Error: Iterated IPC regression paramters have not been estimated.")
      }
    }
  }

  # Select target parameters
  if (is.null(parameter)) {
    parameter <- names(object$InitialIPCRegression)
  }

  # Number of target parameters
  q <- length(parameter)
  param_position <- which(names(object$InitialIPCRegression) %in% parameter)
  significant <- FALSE

  # Add output for initial IPC regression estimates
  if (method_checked == "initial" | method_checked == "both") {

    output[length(output) + 1] <- "Initial IPC regression parameter estimates:"

    for (i in 1:q) {
      output[length(output) + 1] <- paste("Parameter:", parameter[i])
      lm_summary <- capture.output(summary(object$InitialIPCRegression[[param_position[i]]]))
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
  if (method_checked == "iterated" | method_checked == "both") {

    significant <- FALSE

    if(method == "both") {
      output[length(output) + 1] <- ""
    }

    output[length(output) + 1] <- "Iterated IPC regression parameter estimates:"

    for (i in 1:q) {
      output[length(output) + 1] <- paste("Parameter:", parameter[i])
      lm_summary <- capture.output(summary(object$IteratedIPCRegression[[param_position[i]]]))
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
