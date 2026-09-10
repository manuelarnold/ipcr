#' @title Individual Parameter Contribution Regression Summary
#'
#' @description Summary method for objects of class \code{"ipcr"}.
#'
#' @param object An object of class \code{"ipcr"}.
#' @param what A single character string specifying what to display. Options are:
#'   \itemize{
#'     \item \code{"everything"} (default): Displays all summary components.
#'     \item \code{"marginal"}: Shows only marginal linear effects.
#'     \item \code{"ipcr"}: Shows individual parameter contribution regression
#'       coefficients and, when log-variance parameterization was used,
#'       supplementary back-transformed variance effects.
#'   }
#' @param digits A single nonnegative whole number indicating the number of
#'   decimal places to use in the output. Default is 3.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export

summary.ipcr <- function(object, what = "everything", digits = 3, ...) {
  valid_what <- c("everything", "marginal", "ipcr")
  if (!is.character(what) || length(what) != 1L || is.na(what) ||
      !(what %in% valid_what)) {
    stop(
      "'what' must be one of \"everything\", \"marginal\", or \"ipcr\".",
      call. = FALSE
    )
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("'digits' must be a single nonnegative whole number.", call. = FALSE)
  }

  cat("Individual Parameter Contribution Regression\n\n")

  if (inherits(object, "ipcr_it") &&
      identical(object$convergence$likelihood_stable, FALSE)) {
    cat(
      "Likelihood stability warning:\n  ",
      likelihood_stability_warning_ipcr(
        C_l_rate_adjusted = object$convergence$C_l_rate_adjusted,
        threshold = object$convergence$likelihood_warning_threshold
      ),
      "\n\n",
      sep = ""
    )
  }

  if (what %in% c("everything", "marginal")) {
    cat("Marginal Linear Effects:\n\n")

    cat("  MANCOVA:\n")
    if (is.null(object$output$mancova)) {
      cat("    Not available for this result.\n")
    } else {
      # print.summary.Anova.mlm accepts only 1 through 22 for its significant-
      # digits argument. The IPCR tables still use the requested number of
      # decimal places, including zero or values above 22.
      mancova_digits <- max(1, min(digits, 22))
      mancova_output <- utils::capture.output(
        print(
          object$output$mancova,
          digits = mancova_digits,
          SSP = FALSE,
          SSPE = FALSE
        )
      )
      cat(paste0("  ", mancova_output), sep = "\n")
    }
    cat("\n")

    cat("  F-Tests:\n")

    F_tests <- object$output$F_tests
    F_tests_output <- F_tests[, c(
      "Parameter", "Rsquared", "F_Value", "p_Value"
    ), drop = FALSE]
    F_tests_output$stars <- vapply(
      F_tests_output$p_Value,
      add_stars,
      character(1)
    )
    numeric_columns <- c("Rsquared", "F_Value", "p_Value")
    F_tests_output[numeric_columns] <- lapply(
      F_tests_output[numeric_columns],
      round,
      digits = digits
    )
    # Put the headings in the table so that character and numeric columns use
    # the same alignment when the table is printed without column names.
    F_tests_output <- rbind(
      c("Parameter", "Rsquared", "F_Value", "p_Value", ""),
      F_tests_output
    )
    colnames(F_tests_output) <- NULL
    F_tests_output[, 1] <- format(F_tests_output[, 1], justify = "left")
    F_tests_output <- utils::capture.output(
      print(F_tests_output, row.names = FALSE)
    )
    cat(paste0("  ", F_tests_output), sep = "\n")
    cat(paste(
      "   Numerator degrees of freedom (df1):",
      F_tests$DF1[1L]
    ))
    cat("\n")
    cat(paste(
      "   Denominator degrees of freedom (df2):",
      F_tests$DF2[1L]
    ))
    cat("\n\n")
  }

  if (what %in% c("everything", "ipcr")) {
    cat("Individual Parameter Contribution Regression Coefficients:\n")
    ipcr_output <- object$output$ipcr[, c(
      "Parameter", "Predictor", "Estimate", "Std_Error", "t_Value", "p_Value"
    ), drop = FALSE]
    ipcr_output$stars <- vapply(
      ipcr_output$p_Value,
      add_stars,
      character(1)
    )
    numeric_columns <- c("Estimate", "Std_Error", "t_Value", "p_Value")
    ipcr_output[numeric_columns] <- lapply(
      ipcr_output[numeric_columns],
      round,
      digits = digits
    )
    ipcr_output <- rbind(
      c(
        "Parameter", "Predictor", "Estimate", "Std_Error", "t_Value",
        "p_Value", ""
      ),
      ipcr_output
    )
    colnames(ipcr_output) <- NULL
    ipcr_output[, c(1, 2)] <- lapply(
      ipcr_output[, c(1, 2), drop = FALSE],
      format,
      justify = "left"
    )
    ipcr_output <- utils::capture.output(
      print(ipcr_output, row.names = FALSE)
    )
    cat(paste0("  ", ipcr_output), sep = "\n")

    natural_scale_effects <- object$output$natural_scale_effects
    if (!is.null(natural_scale_effects) &&
        NROW(natural_scale_effects) > 0L) {
      cat("\n\nBack-Transformed Variance Effects:\n")
      natural_output <- natural_scale_effects[, c(
        "Parameter", "Predictor", "Log_Effect", "Variance_Ratio",
        "Reference_Variance", "Variance_Change"
      ), drop = FALSE]
      numeric_columns <- c(
        "Log_Effect", "Variance_Ratio", "Reference_Variance",
        "Variance_Change"
      )
      natural_output[numeric_columns] <- lapply(
        natural_output[numeric_columns],
        round,
        digits = digits
      )
      natural_output <- utils::capture.output(
        print(natural_output, row.names = FALSE)
      )
      cat(paste0("  ", natural_output), sep = "\n")
      cat(
        "\n  Variance changes compare a design-matrix column value of 0 ",
        "with 1 while holding all other non-intercept columns at 0."
      )
    }
  }

  invisible(NULL)
}
