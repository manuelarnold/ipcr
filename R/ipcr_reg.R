#' Regularized Individual Parameter Contribution Regression
#'
#' Applies a cross-validated multivariate group elastic net to the IPCs and
#' predictor design matrix stored in an existing [ipcr()] or [ipcr_it()]
#' result. Regularization is confined to the second-stage IPC regression;
#' scores and IPCs are not recomputed.
#'
#' @param x An object inheriting from class `"ipcr"`, including an
#'   `"ipcr_it"` result.
#' @param alpha A single number in `[0, 1]` controlling the elastic-net
#'   mixture. The default, `0.5`, balances the sparsity of the group lasso with
#'   ridge stabilization for correlated predictors. `alpha = 1` is group
#'   lasso and `alpha = 0` is group ridge regression.
#' @param s A single selection rule or penalty value passed to the coefficient
#'   and prediction methods for the cross-validated `glmnet` fit. Supported
#'   character values are `"lambda.1se"` (default) and `"lambda.min"`; a
#'   positive numeric lambda is also accepted.
#' @param nfolds A whole number of cross-validation folds. The default is 10.
#'   If fewer than 10 complete IPC rows are available, the number is reduced
#'   to the number of available rows with a warning.
#' @param foldid Optional vector assigning each complete IPC row to a
#'   cross-validation fold. At least three folds are required. If omitted, a
#'   balanced random assignment is generated and stored in the result.
#' @param standardize Logical; whether `glmnet` should standardize predictor
#'   design columns. The default is `TRUE`.
#' @param standardize_response Logical; whether multiple IPC response columns
#'   should be standardized before fitting. The default is `TRUE`, which keeps
#'   parameters measured on different scales from dominating the joint loss.
#'   This argument is ignored when there is only one IPC parameter.
#' @param ... Additional named arguments passed to [glmnet::cv.glmnet()]. The
#'   arguments controlled explicitly by `ipcr_reg()` cannot be supplied here.
#'
#' @details
#' With two or more IPC parameters, `ipcr_reg()` uses
#' `family = "mgaussian"`. The coefficient vector for each predictor design
#' column is penalized jointly across all IPC responses, so that column is
#' selected for or excluded from the second-stage IPC regression as a whole.
#' With one IPC parameter, the function uses `family = "gaussian"`. At least
#' two non-intercept predictor design columns are required because `glmnet`
#' does not fit a one-column design matrix.
#'
#' Factors are represented by the same design-matrix contrasts used by the
#' original IPC regression. Each contrast is penalized separately. The result
#' records both contrast-level selections and their originating predictor
#' terms; consequently, selection for a factor with more than two levels can
#' depend on its contrast coding and reference level.
#'
#' The selected coefficients are intended for prediction and exploratory
#' variable selection. Ordinary standard errors, confidence intervals, and
#' p-values are not available for the regularized fit. Cross-validation
#' assesses the second-stage prediction of the already-computed IPCs; it does
#' not refit or cross-validate the first-stage statistical model.
#'
#' `ipcr_it()` returns an ordinary result only after convergence. For safety,
#' `ipcr_reg()` also rejects legacy or manually constructed iterated results
#' whose convergence status is not `TRUE`.
#'
#' @section Methods:
#' `coef()` and `coefficients()` return the regularized coefficient matrix.
#' `predict()` returns predicted IPCs, optionally for `newdata`; `fitted()`
#' returns predictions for the analyzed rows, `residuals()` returns the
#' corresponding IPC residuals, and `nobs()` returns the number of analyzed
#' rows. Because `glmnet` does not provide ordinary second-stage sampling
#' covariance estimates, `confint()`, `vcov()`, `effects()`, and `sigma()` are
#' deliberately unavailable for regularized IPCR results.
#'
#' @return An object whose class begins with `"ipcr_reg"` and which otherwise
#'   retains the input object's classes and components. The additional
#'   `regularization` component contains the `cv.glmnet` fit, selected lambda,
#'   coefficient matrix, selected predictor contrasts and terms, standardized
#'   predictor importance, fold assignments, and design-column metadata.
#'
#' @examples
#' set.seed(1)
#' n <- 80L
#' predictors <- data.frame(
#'   group = rep(0:1, each = n / 2),
#'   age = rnorm(n),
#'   noise = rnorm(n)
#' )
#' x <- rnorm(n)
#' y <- 1 + x + 0.4 * predictors$group * x + rnorm(n)
#' fit <- lm(y ~ x)
#' standard <- ipcr(fit, predictors)
#' regularized <- ipcr_reg(standard, nfolds = 5)
#' summary(regularized)
#'
#' @seealso [ipcr()], [ipcr_it()], [glmnet::cv.glmnet()]
#' @export
ipcr_reg <- function(
    x, alpha = 0.5, s = "lambda.1se", nfolds = 10L, foldid = NULL,
    standardize = TRUE, standardize_response = TRUE, ...) {
  components <- validate_ipcr_reg_object(x)
  validate_ipcr_reg_scalar_arguments(
    alpha = alpha,
    s = s,
    nfolds = nfolds,
    standardize = standardize,
    standardize_response = standardize_response
  )

  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names) || any(is.na(dot_names) | !nzchar(dot_names))) {
      stop("Additional arguments in '...' must be named.", call. = FALSE)
    }
    protected <- c(
      "x", "y", "family", "alpha", "nfolds", "foldid", "standardize",
      "standardize.response", "intercept"
    )
    duplicated_arguments <- intersect(dot_names, protected)
    if (length(duplicated_arguments) > 0L) {
      stop(
        "The following arguments are controlled by ipcr_reg() and cannot ",
        "be supplied through '...': ",
        paste(duplicated_arguments, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (inherits(x, "ipcr_it") &&
      !isTRUE(x$convergence$converged)) {
    stop(
      "Regularized IPCR requires a converged ipcr_it object.",
      call. = FALSE
    )
  }
  if (alpha == 0) {
    warning(
      "alpha = 0 fits group ridge regression and does not generally select ",
      "predictors by setting their coefficients to zero.",
      call. = FALSE
    )
  }

  fold_components <- prepare_ipcr_reg_folds(
    n = NROW(components$responses),
    nfolds = nfolds,
    foldid = foldid
  )
  family <- if (NCOL(components$responses) > 1L) {
    "mgaussian"
  } else {
    "gaussian"
  }
  response <- if (identical(family, "mgaussian")) {
    components$responses
  } else {
    drop(components$responses)
  }

  fit_arguments <- list(
    x = components$design,
    y = response,
    family = family,
    alpha = alpha,
    nfolds = fold_components$nfolds,
    foldid = fold_components$foldid,
    standardize = standardize,
    intercept = TRUE
  )
  if (identical(family, "mgaussian")) {
    fit_arguments$standardize.response <- standardize_response
  }
  fit_arguments <- c(fit_arguments, dots)

  regularized_fit <- tryCatch(
    do.call(glmnet::cv.glmnet, fit_arguments),
    error = function(e) {
      stop(
        "The regularized IPC regression could not be fitted: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  coefficient_matrix <- extract_ipcr_reg_coefficients(
    fit = regularized_fit,
    s = s,
    parameters = colnames(components$responses),
    design_columns = colnames(components$design),
    family = family
  )
  output <- regularized_output_ipcr(
    coefficients = coefficient_matrix,
    design = components$design,
    responses = components$responses,
    design_terms = components$design_terms
  )
  selected_lambda <- resolve_ipcr_reg_lambda(regularized_fit, s)

  result <- x
  result$info$regularized <- TRUE
  result$regularization <- list(
    fit = regularized_fit,
    family = family,
    alpha = alpha,
    s = s,
    lambda = selected_lambda,
    nfolds = fold_components$nfolds,
    foldid = fold_components$foldid,
    standardize = standardize,
    standardize_response = if (identical(family, "mgaussian")) {
      standardize_response
    } else {
      NULL
    },
    coefficients = coefficient_matrix,
    selected_columns = output$selected_columns,
    selected_predictors = output$selected_predictors,
    importance = output$importance,
    predictor_terms = output$predictor_terms,
    coefficient_table = output$coefficient_table,
    design_columns = colnames(components$design),
    design_terms = components$design_terms,
    call = match.call()
  )
  result$output$regularized_coefficients <- output$coefficient_table
  result$output$regularized_importance <- output$importance
  result$output$regularized_predictors <- output$predictor_terms
  class(result) <- unique(c("ipcr_reg", class(x)))
  result
}

validate_ipcr_reg_object <- function(x) {
  if (!inherits(x, "ipcr")) {
    stop(
      "'x' must be an object returned by ipcr() or ipcr_it().",
      call. = FALSE
    )
  }
  if (inherits(x, "ipcr_reg")) {
    stop(
      "'x' is already regularized; supply its original ipcr or ipcr_it object.",
      call. = FALSE
    )
  }
  if (!is.matrix(x$IPCs) || !is.numeric(x$IPCs) ||
      NROW(x$IPCs) < 1L || NCOL(x$IPCs) < 1L) {
    stop("The IPCR object does not contain a usable numeric IPC matrix.", call. = FALSE)
  }
  if (is.null(x$complete_rows) || !is.logical(x$complete_rows) ||
      length(x$complete_rows) != NROW(x$IPCs) || anyNA(x$complete_rows)) {
    stop(
      "The IPCR object does not contain a valid 'complete_rows' indicator.",
      call. = FALSE
    )
  }
  if (!any(x$complete_rows)) {
    stop("The IPCR object contains no complete second-stage rows.", call. = FALSE)
  }
  if (!inherits(x$mlm, "lm")) {
    stop("The IPCR object does not contain its second-stage regression.", call. = FALSE)
  }

  responses <- x$IPCs[x$complete_rows, , drop = FALSE]
  parameter_names <- x$info$parameters
  if (is.null(parameter_names) || length(parameter_names) != NCOL(responses)) {
    stop("Unable to align IPC responses with the stored parameter names.", call. = FALSE)
  }
  colnames(responses) <- parameter_names
  if (any(!is.finite(responses))) {
    stop("The complete-case IPC responses must be finite.", call. = FALSE)
  }
  response_scales <- apply(responses, 2L, function(column) {
    sqrt(mean((column - mean(column))^2))
  })
  if (any(!is.finite(response_scales) | response_scales <= 0)) {
    invalid_response <- !is.finite(response_scales) | response_scales <= 0
    stop(
      "Every IPC response must have positive variation for regularization: ",
      paste(parameter_names[invalid_response], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  full_design <- stats::model.matrix(x$mlm)
  design_assign <- attr(full_design, "assign")
  terms <- stats::terms(x$mlm)
  term_labels <- attr(terms, "term.labels")
  predictor_columns <- which(design_assign != 0L)
  design <- full_design[, predictor_columns, drop = FALSE]
  if (NCOL(design) < 2L) {
    stop(
      "Regularized IPCR currently requires at least two predictor design ",
      "columns because glmnet does not fit a one-column design matrix.",
      call. = FALSE
    )
  }
  if (NROW(design) != NROW(responses)) {
    stop(
      "Unable to align the stored IPC responses with the second-stage design matrix.",
      call. = FALSE
    )
  }
  if (!is.numeric(design) || any(!is.finite(design))) {
    stop("The second-stage predictor design matrix must be finite and numeric.", call. = FALSE)
  }
  if (is.null(design_assign) || any(design_assign[predictor_columns] < 1L) ||
      any(design_assign[predictor_columns] > length(term_labels))) {
    stop("Unable to map predictor design columns to their model terms.", call. = FALSE)
  }
  design_terms <- unname(term_labels[design_assign[predictor_columns]])
  names(design_terms) <- colnames(design)

  model_response <- as.matrix(stats::model.response(stats::model.frame(x$mlm)))
  if (!identical(dim(model_response), dim(responses)) ||
      !isTRUE(all.equal(
        unname(model_response), unname(responses),
        tolerance = sqrt(.Machine$double.eps)
      ))) {
    stop(
      "The stored IPC matrix is inconsistent with the second-stage regression response.",
      call. = FALSE
    )
  }

  list(
    responses = responses,
    design = design,
    design_terms = design_terms
  )
}

validate_ipcr_reg_scalar_arguments <- function(
    alpha, s, nfolds, standardize, standardize_response) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single finite number between zero and one.", call. = FALSE)
  }
  validate_ipcr_reg_s(s)
  if (!is.numeric(nfolds) || length(nfolds) != 1L || is.na(nfolds) ||
      !is.finite(nfolds) || nfolds < 3 || nfolds != floor(nfolds)) {
    stop("'nfolds' must be a single whole number of at least three.", call. = FALSE)
  }
  for (argument in c("standardize", "standardize_response")) {
    value <- get(argument)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("'", argument, "' must be either TRUE or FALSE.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

validate_ipcr_reg_s <- function(s) {
  valid_character <- is.character(s) && length(s) == 1L && !is.na(s) &&
    s %in% c("lambda.1se", "lambda.min")
  valid_numeric <- is.numeric(s) && length(s) == 1L && !is.na(s) &&
    is.finite(s) && s > 0
  if (!valid_character && !valid_numeric) {
    stop(
      "'s' must be \"lambda.1se\", \"lambda.min\", or a single positive lambda.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

prepare_ipcr_reg_folds <- function(n, nfolds, foldid) {
  if (n < 3L) {
    stop("Regularized IPCR requires at least three complete IPC rows.", call. = FALSE)
  }
  if (is.null(foldid)) {
    if (nfolds > n) {
      warning(
        "'nfolds' exceeds the number of complete IPC rows and was reduced from ",
        nfolds, " to ", n, ".",
        call. = FALSE
      )
      nfolds <- n
    }
    nfolds <- as.integer(nfolds)
    foldid <- sample(rep(seq_len(nfolds), length.out = n))
  } else {
    if (!is.numeric(foldid) || length(foldid) != n || anyNA(foldid) ||
        any(!is.finite(foldid)) || any(foldid < 1) ||
        any(foldid != floor(foldid))) {
      stop(
        "'foldid' must contain one positive whole-number fold assignment for ",
        "each complete IPC row.",
        call. = FALSE
      )
    }
    foldid <- match(foldid, unique(foldid))
    nfolds <- length(unique(foldid))
    if (nfolds < 3L) {
      stop("'foldid' must define at least three folds.", call. = FALSE)
    }
  }
  list(nfolds = as.integer(nfolds), foldid = as.integer(foldid))
}

extract_ipcr_reg_coefficients <- function(
    fit, s, parameters, design_columns, family) {
  validate_ipcr_reg_s(s)
  raw_coefficients <- stats::coef(fit, s = s)
  expected_rows <- c("(Intercept)", design_columns)

  if (identical(family, "mgaussian")) {
    if (!is.list(raw_coefficients) || length(raw_coefficients) != length(parameters)) {
      stop("Unable to extract the multivariate regularized coefficients.", call. = FALSE)
    }
    coefficient_matrix <- do.call(cbind, lapply(raw_coefficients, function(value) {
      as.numeric(value[, 1L])
    }))
    rownames(coefficient_matrix) <- rownames(raw_coefficients[[1L]])
    colnames(coefficient_matrix) <- parameters
  } else {
    coefficient_matrix <- matrix(
      as.numeric(raw_coefficients[, 1L]),
      ncol = 1L,
      dimnames = list(rownames(raw_coefficients), parameters)
    )
  }
  if (!identical(rownames(coefficient_matrix), expected_rows)) {
    stop(
      "Unable to align regularized coefficients with predictor design columns.",
      call. = FALSE
    )
  }
  coefficient_matrix
}

resolve_ipcr_reg_lambda <- function(fit, s) {
  if (is.character(s)) {
    return(as.numeric(fit[[s]]))
  }
  as.numeric(s)
}

regularized_output_ipcr <- function(
    coefficients, design, responses, design_terms) {
  design_columns <- colnames(design)
  predictor_coefficients <- coefficients[design_columns, , drop = FALSE]
  selected <- apply(predictor_coefficients != 0, 1L, any)

  x_scales <- apply(design, 2L, function(column) {
    sqrt(mean((column - mean(column))^2))
  })
  y_scales <- apply(responses, 2L, function(column) {
    sqrt(mean((column - mean(column))^2))
  })
  standardized_coefficients <- sweep(
    predictor_coefficients,
    1L,
    x_scales,
    FUN = "*"
  )
  standardized_coefficients <- sweep(
    standardized_coefficients,
    2L,
    y_scales,
    FUN = "/"
  )
  standardized_coefficients[!is.finite(standardized_coefficients)] <- 0
  column_importance <- sqrt(rowSums(standardized_coefficients^2))

  importance <- data.frame(
    Design_Column = design_columns,
    Predictor = unname(design_terms),
    Importance = unname(column_importance),
    Selected = unname(selected),
    stringsAsFactors = FALSE
  )
  term_indices <- split(seq_along(design_terms), factor(
    design_terms,
    levels = unique(design_terms)
  ))
  predictor_terms <- do.call(rbind, lapply(names(term_indices), function(term) {
    indices <- term_indices[[term]]
    selected_columns <- design_columns[indices][selected[indices]]
    data.frame(
      Predictor = term,
      Design_Columns = paste(design_columns[indices], collapse = ", "),
      Selected_Columns = paste(selected_columns, collapse = ", "),
      Importance = sqrt(sum(column_importance[indices]^2)),
      Selected = any(selected[indices]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(predictor_terms) <- NULL

  coefficient_grid <- expand.grid(
    Design_Column = rownames(coefficients),
    Parameter = colnames(coefficients),
    stringsAsFactors = FALSE
  )
  coefficient_grid$Predictor <- c(
    "(Intercept)", unname(design_terms)
  )[match(coefficient_grid$Design_Column, rownames(coefficients))]
  coefficient_grid$Estimate <- as.vector(coefficients)
  coefficient_grid$Selected <- ifelse(
    coefficient_grid$Design_Column == "(Intercept)",
    NA,
    coefficient_grid$Estimate != 0
  )
  coefficient_table <- coefficient_grid[, c(
    "Parameter", "Design_Column", "Predictor", "Estimate", "Selected"
  )]

  list(
    coefficient_table = coefficient_table,
    importance = importance,
    predictor_terms = predictor_terms,
    selected_columns = design_columns[selected],
    selected_predictors = predictor_terms$Predictor[predictor_terms$Selected]
  )
}

regularized_design_matrix_ipcr <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    full_design <- stats::model.matrix(object$mlm)
  } else {
    predictor_terms <- stats::delete.response(stats::terms(object$mlm))
    full_design <- stats::model.matrix(
      predictor_terms,
      data = newdata,
      contrasts.arg = object$mlm$contrasts,
      xlev = object$mlm$xlevels
    )
  }
  required <- object$regularization$design_columns
  if (!all(required %in% colnames(full_design))) {
    stop(
      "The prediction data do not produce all required predictor design columns.",
      call. = FALSE
    )
  }
  full_design[, required, drop = FALSE]
}

#' @exportS3Method stats::coef
coef.ipcr_reg <- function(object, s = object$regularization$s, ...) {
  extract_ipcr_reg_coefficients(
    fit = object$regularization$fit,
    s = s,
    parameters = object$info$parameters,
    design_columns = object$regularization$design_columns,
    family = object$regularization$family
  )
}

#' @export
coefficients.ipcr_reg <- function(object, s = object$regularization$s, ...) {
  coef.ipcr_reg(object, s = s, ...)
}

#' @export
predict.ipcr_reg <- function(
    object, newdata = NULL, s = object$regularization$s, ...) {
  validate_ipcr_reg_s(s)
  dots <- list(...)
  if ("type" %in% names(dots)) {
    stop("'type' is fixed to \"response\" for ipcr_reg predictions.", call. = FALSE)
  }
  design <- regularized_design_matrix_ipcr(object, newdata = newdata)
  prediction <- do.call(
    stats::predict,
    c(list(
      object = object$regularization$fit,
      newx = design,
      s = s,
      type = "response"
    ), dots)
  )
  if (identical(object$regularization$family, "mgaussian")) {
    prediction <- prediction[, , 1L, drop = FALSE]
    prediction <- matrix(
      prediction,
      nrow = NROW(design),
      ncol = length(object$info$parameters),
      dimnames = list(rownames(design), object$info$parameters)
    )
  } else {
    prediction <- matrix(
      as.numeric(prediction[, 1L]),
      ncol = 1L,
      dimnames = list(rownames(design), object$info$parameters)
    )
  }
  prediction
}

#' @exportS3Method stats::fitted
fitted.ipcr_reg <- function(object, ...) {
  predict.ipcr_reg(object, ...)
}

#' @export
residuals.ipcr_reg <- function(object, ...) {
  responses <- object$IPCs[object$complete_rows, , drop = FALSE]
  responses - fitted.ipcr_reg(object, ...)
}

#' @export
nobs.ipcr_reg <- function(object, ...) {
  as.integer(sum(object$complete_rows))
}

regularized_inference_unavailable_ipcr <- function(generic) {
  stop(
    generic,
    " is not available for regularized IPCR because glmnet does not provide ",
    "ordinary second-stage sampling covariance estimates.",
    call. = FALSE
  )
}

#' @export
confint.ipcr_reg <- function(object, ...) {
  regularized_inference_unavailable_ipcr("confint()")
}

#' @export
vcov.ipcr_reg <- function(object, ...) {
  regularized_inference_unavailable_ipcr("vcov()")
}

#' @export
effects.ipcr_reg <- function(object, ...) {
  regularized_inference_unavailable_ipcr("effects()")
}

#' @export
sigma.ipcr_reg <- function(object, ...) {
  regularized_inference_unavailable_ipcr("sigma()")
}

#' @export
print.ipcr_reg <- function(x, ...) {
  if (inherits(x, "ipcr_it")) {
    print.ipcr_it(x, ...)
  } else {
    print.ipcr(x, ...)
  }
  method <- if (identical(x$regularization$family, "mgaussian")) {
    "Multivariate group elastic net"
  } else {
    "Elastic net"
  }
  cat("Regularization:", method, "\n")
  cat("Elastic-net alpha:", x$regularization$alpha, "\n")
  cat("Selected lambda: ", format(x$regularization$lambda), " (",
      x$regularization$s, ")\n", sep = "")
  selected <- x$regularization$selected_predictors
  cat(
    "Selected predictors:",
    if (length(selected) > 0L) paste(selected, collapse = ", ") else "none",
    "\n"
  )
  invisible(x)
}

#' Summarize a Regularized IPCR Result
#'
#' @param object An object returned by [ipcr_reg()].
#' @param what A single character string: `"everything"` (default),
#'   `"predictors"`, or `"coefficients"`.
#' @param digits A nonnegative whole number of digits used for printing.
#' @param ... Additional arguments passed to or from methods.
#'
#' @return Invisibly returns `NULL`.
#' @export
summary.ipcr_reg <- function(
    object, what = "everything", digits = 3, ...) {
  valid_what <- c("everything", "predictors", "coefficients")
  if (!is.character(what) || length(what) != 1L || is.na(what) ||
      !(what %in% valid_what)) {
    stop(
      "'what' must be one of \"everything\", \"predictors\", or \"coefficients\".",
      call. = FALSE
    )
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("'digits' must be a single nonnegative whole number.", call. = FALSE)
  }

  cat("Regularized Individual Parameter Contribution Regression\n\n")

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

  cat("Method:", if (identical(object$regularization$family, "mgaussian")) {
    "Multivariate group elastic net"
  } else {
    "Elastic net"
  }, "\n")
  cat("Alpha:", object$regularization$alpha, "\n")
  cat("Lambda: ", round(object$regularization$lambda, digits),
      " (", object$regularization$s, ")\n", sep = "")
  cat("Cross-validation folds:", object$regularization$nfolds, "\n\n")

  if (what %in% c("everything", "predictors")) {
    cat("Predictor selection and standardized importance:\n")
    predictor_table <- object$regularization$predictor_terms
    predictor_table$Importance <- round(predictor_table$Importance, digits)
    print(predictor_table, row.names = FALSE)
    cat("\n")
  }
  if (what %in% c("everything", "coefficients")) {
    cat("Regularized IPC regression coefficients:\n")
    coefficient_table <- object$regularization$coefficient_table
    coefficient_table$Estimate <- round(coefficient_table$Estimate, digits)
    print(coefficient_table, row.names = FALSE)
    cat("\nOrdinary standard errors and p-values are not available.\n")
  }
  invisible(NULL)
}
