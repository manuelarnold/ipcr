check_arguments_ipcr <- function(x, predictors, analytic) {

  # Check fitted model ----
  check_model_support_ipcr(x)
  check_model_scope_ipcr(x)

  # Check argument analytic ----

  ## Check if analytic argument
  if (!is.logical(analytic) || length(analytic) != 1L || is.na(analytic)) {
    stop("'analytic' must be either TRUE or FALSE.", call. = FALSE)
  }

  normalize_analytic_ipcr(x, analytic)
}

validate_ipcr_it_common_arguments <- function(
    analytic, conv, learning_rate, max_it, iteration_info) {
  if (!is.logical(analytic) || length(analytic) != 1L || is.na(analytic)) {
    stop("'analytic' must be either TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conv) || length(conv) != 1L || is.na(conv) ||
      !is.finite(conv) || conv <= 0) {
    stop("'conv' must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(learning_rate) || length(learning_rate) == 0L ||
      anyNA(learning_rate) || any(!is.finite(learning_rate)) ||
      any(learning_rate <= 0 | learning_rate > 1)) {
    stop(
      "'learning_rate' must contain finite values greater than 0 and no greater than 1.",
      call. = FALSE
    )
  }
  if (!is.numeric(max_it) || length(max_it) != 1L || is.na(max_it) ||
      !is.finite(max_it) || max_it < 1 || max_it != floor(max_it)) {
    stop("'max_it' must be one positive whole number.", call. = FALSE)
  }
  if (!is.logical(iteration_info) || length(iteration_info) != 1L ||
      is.na(iteration_info)) {
    stop("'iteration_info' must be either TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}

prepare_predictors_ipcr <- function(predictors, n) {
  vector_input <- is.atomic(predictors) && is.null(dim(predictors))
  supported_container <- vector_input || is.matrix(predictors) ||
    is.data.frame(predictors)

  if (!supported_container) {
    stop(
      "Argument 'predictors' must be a vector, matrix, or data.frame.",
      call. = FALSE
    )
  }

  predictors <- tryCatch(
    as.data.frame(predictors, optional = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(
        "Argument 'predictors' could not be converted to a data.frame: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (NCOL(predictors) == 0L) {
    stop(
      "Argument 'predictors' must contain at least one column.",
      call. = FALSE
    )
  }
  if (NROW(predictors) != n) {
    stop(
      "The number of rows in 'predictors' (", NROW(predictors),
      ") must equal the number of IPC units in the fitted model (", n,
      "). Each predictor row must correspond, in order, to the same IPC unit ",
      "used to fit the model.",
      call. = FALSE
    )
  }

  predictor_names <- names(predictors)
  if (is.null(predictor_names)) {
    predictor_names <- rep("", NCOL(predictors))
  }
  unnamed <- is.na(predictor_names) | !nzchar(predictor_names)
  if (anyDuplicated(predictor_names[!unnamed])) {
    stop("Predictor names must be unique.", call. = FALSE)
  }
  if (any(unnamed)) {
    used_names <- predictor_names[!unnamed]
    next_number <- 1L
    for (column in which(unnamed)) {
      candidate <- paste0("predictor", next_number)
      while (candidate %in% used_names) {
        next_number <- next_number + 1L
        candidate <- paste0("predictor", next_number)
      }
      predictor_names[column] <- candidate
      used_names <- c(used_names, candidate)
      next_number <- next_number + 1L
    }
    names(predictors) <- predictor_names
    warning(
      "Unnamed predictor columns were assigned the names: ",
      paste(predictor_names[unnamed], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  supported_column <- vapply(
    predictors,
    function(column) {
      is.numeric(column) || is.logical(column) || is.factor(column)
    },
    logical(1)
  )
  if (any(!supported_column)) {
    stop(
      "Unsupported predictor columns: ",
      paste(names(predictors)[!supported_column], collapse = ", "),
      ". Predictors must be numeric, logical, or factors; convert character ",
      "columns to factors or numeric indicators.",
      call. = FALSE
    )
  }

  finite_numeric <- vapply(
    predictors,
    function(column) {
      if (!is.numeric(column)) {
        return(TRUE)
      }
      all(is.finite(column[!is.na(column)]))
    },
    logical(1)
  )
  if (any(!finite_numeric)) {
    stop(
      "Numeric predictor columns must not contain infinite values: ",
      paste(names(predictors)[!finite_numeric], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  complete_rows <- stats::complete.cases(predictors)
  if (!any(complete_rows)) {
    stop(
      "Argument 'predictors' has no complete rows available for regression.",
      call. = FALSE
    )
  }
  if (sum(complete_rows) < 2L) {
    stop(
      "At least two complete predictor rows are required for regression.",
      call. = FALSE
    )
  }
  if (any(!complete_rows)) {
    warning(
      sum(!complete_rows),
      " row", if (sum(!complete_rows) == 1L) "" else "s",
      " with missing predictor values will be omitted from the IPC regression; ",
      "the corresponding IPCs and scores remain in the returned object.",
      call. = FALSE
    )
  }

  analysis_predictors <- droplevels(
    predictors[complete_rows, , drop = FALSE]
  )
  factor_levels <- vapply(
    analysis_predictors,
    function(column) !is.factor(column) || nlevels(column) >= 2L,
    logical(1)
  )
  if (any(!factor_levels)) {
    stop(
      "Factor predictors must have at least two observed levels after rows ",
      "with missing values are omitted: ",
      paste(names(analysis_predictors)[!factor_levels], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    predictors = predictors,
    analysis_predictors = analysis_predictors,
    complete_rows = complete_rows
  )
}
