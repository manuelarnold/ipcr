# Classify fitted models without preventing generic sandwich-compatible
# classes from being attempted.
model_support_ipcr <- function(x) {
  if (inherits(x, "lavaan") ||
      inherits(x, "MxModel") ||
      inherits(x, "lmerMod") ||
      inherits(x, "glmerMod") ||
      inherits(x, "glm") ||
      identical(class(x)[1], "lm")) {
    return("official")
  }

  "experimental"
}

check_model_support_ipcr <- function(x) {
  support <- model_support_ipcr(x)
  if (identical(support, "experimental")) {
    warning(
      "Models of class '",
      paste(class(x), collapse = "/"),
      "' are not officially supported. ipcr will attempt the generic ",
      "sandwich score and bread methods; validate the parameterization ",
      "carefully.",
      call. = FALSE
    )
  }
  invisible(support)
}

mixed_model_group_ipcr <- function(x) {
  grouping_factors <- lme4::getME(x, "flist")
  if (length(grouping_factors) != 1L) {
    stop(
      "Official mixed-model support currently requires exactly one grouping factor.",
      call. = FALSE
    )
  }
  droplevels(grouping_factors[[1]])
}

check_model_scope_ipcr <- function(x) {
  if (inherits(x, "lavaan")) {
    check_lavaan_scope_ipcr(x)
  }

  if (inherits(x, "lmerMod")) {
    mixed_model_group_ipcr(x)
    if (lme4::isREML(x)) {
      stop(
        "Official lmer support currently requires a model fitted with REML = FALSE.",
        call. = FALSE
      )
    }
  }

  if (inherits(x, "glmerMod")) {
    mixed_model_group_ipcr(x)
    model_family <- stats::family(x)$family
    if (!model_family %in% c("binomial", "poisson")) {
      stop(
        "Official glmer support is currently limited to binomial and Poisson models.",
        call. = FALSE
      )
    }
  }

  if (inherits(x, "merMod")) {
    if (!is.null(x@call$weights)) {
      stop(
        "Official mixed-model support does not currently include prior weights.",
        call. = FALSE
      )
    }
    model_formula <- paste(deparse(stats::formula(x)), collapse = " ")
    if (grepl("cbind\\s*\\(", model_formula)) {
      stop(
        "Official mixed-model support does not currently include cbind responses.",
        call. = FALSE
      )
    }
  }

  if (inherits(x, "MxModel")) {
    openmx_raw_data_ipcr(x)
  }

  invisible(TRUE)
}

check_lavaan_scope_ipcr <- function(x) {
  estimator <- toupper(x@Options$estimator)
  estimator_original <- toupper(x@Options$estimator.orig)
  categorical <- isTRUE(x@Model@categorical)
  supported_estimators <- c("ML", "GLS", "ULS", "WLS")

  if (categorical && !identical(estimator, "WLS")) {
    stop(
      "Ordered-indicator lavaan models currently require estimator = 'WLS' ",
      "because estimator '",
      estimator_original,
      "' does not provide both compatible casewise scores and expected ",
      "information.",
      call. = FALSE
    )
  }

  if (!estimator %in% supported_estimators) {
    stop(
      "The lavaan estimator '",
      estimator_original,
      "' does not provide both compatible casewise scores and expected ",
      "information for ipcr.",
      call. = FALSE
    )
  }

  if (!identical(x@Options$missing, "listwise") &&
      !identical(estimator, "ML")) {
    stop(
      "Incomplete-data lavaan models are supported only when lavaan uses ",
      "full-information ML.",
      call. = FALSE
    )
  }

  if (x@Data@nlevels != 1L) {
    stop(
      "Multilevel lavaan models are not currently supported by ipcr.",
      call. = FALSE
    )
  }

  if (isTRUE(x@Model@conditional.x) &&
      !(identical(estimator, "ML") && !categorical)) {
    stop(
      "This conditional.x lavaan model does not provide compatible casewise scores.",
      call. = FALSE
    )
  }

  constraint_rows <- NROW(x@Model@con.jac)
  simple_equality <- isTRUE(x@Model@eq.constraints) ||
    isTRUE(x@Model@ceq.simple.only)
  if (constraint_rows > 0L && !simple_equality) {
    stop(
      "Only simple linear equality constraints are supported for lavaan models; ",
      "nonlinear and inequality constraints are not supported.",
      call. = FALSE
    )
  }
  if (simple_equality) {
    lavaan_constraint_reduction_ipcr(x@Model)
  }

  invisible(TRUE)
}

normalize_analytic_ipcr <- function(x, analytic) {
  if (!inherits(x, "MxModel")) {
    return(analytic)
  }

  if (isTRUE(analytic) && !inherits(x, "MxRAMModel")) {
    warning(
      "Analytical scores are unavailable for non-RAM MxModels; using ",
      "numerical OpenMx row scores instead.",
      call. = FALSE
    )
    return(FALSE)
  }

  if (isTRUE(analytic)) {
    issue <- openmx_analytical_score_issue_ipcr(x)
    if (!is.null(issue)) {
      warning(
        "Analytical OpenMx RAM scores do not support ",
        issue,
        "; using numerical OpenMx row scores instead.",
        call. = FALSE
      )
      return(FALSE)
    }
  }

  analytic
}

prepare_model_ipcr <- function(x, analytic) {
  if (!is.logical(analytic) || length(analytic) != 1L || is.na(analytic)) {
    stop("'analytic' must be either TRUE or FALSE.", call. = FALSE)
  }
  check_model_support_ipcr(x)
  check_model_scope_ipcr(x)
  normalize_analytic_ipcr(x, analytic)
}

openmx_raw_data_ipcr <- function(x) {
  data_type <- tryCatch(x$data$type, error = function(e) NULL)
  observed <- tryCatch(x$data$observed, error = function(e) NULL)
  if (!identical(data_type, "raw") || is.null(observed)) {
    stop(
      "OpenMx support currently requires raw data stored in the fitted model.",
      call. = FALSE
    )
  }
  observed
}

# Validate fitted-model convergence before extracting IPC components.
check_model_convergence_ipcr <- function(x, ...) {
  UseMethod("check_model_convergence_ipcr")
}

#' @noRd
#' @export
check_model_convergence_ipcr.default <- function(x, ...) {
  invisible(TRUE)
}

#' @noRd
#' @export
check_model_convergence_ipcr.glm <- function(x, ...) {
  if (isFALSE(x$converged)) {
    stop("The fitted glm model did not converge.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
#' @export
check_model_convergence_ipcr.nls <- function(x, ...) {
  if (!is.null(x$convInfo$isConv) && isFALSE(x$convInfo$isConv)) {
    stop("The fitted nls model did not converge.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
#' @export
check_model_convergence_ipcr.lavaan <- function(x, ...) {
  converged <- tryCatch(
    lavaan::lavInspect(x, what = "converged"),
    error = function(e) {
      stop(
        "Unable to verify whether the fitted lavaan model converged: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!isTRUE(converged)) {
    stop("The fitted lavaan model did not converge.", call. = FALSE)
  }
  invisible(TRUE)
}

check_mermod_convergence_ipcr <- function(x) {
  optimizer_code <- unlist(x@optinfo$conv$opt, use.names = FALSE)
  convergence_messages <- x@optinfo$conv$lme4$messages

  optimizer_failed <- length(optimizer_code) > 0L &&
    any(is.na(optimizer_code) | optimizer_code != 0)
  messages_present <- length(convergence_messages) > 0L

  if (optimizer_failed || messages_present) {
    details <- if (messages_present) {
      paste(convergence_messages, collapse = "; ")
    } else {
      paste("optimizer code", paste(optimizer_code, collapse = ", "))
    }
    stop(
      "The fitted lme4 model did not converge: ",
      details,
      call. = FALSE
    )
  }

  if (lme4::isSingular(x)) {
    warning(
      "The fitted lme4 model has a singular random-effects structure.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
#' @export
check_model_convergence_ipcr.lmerMod <- function(x, ...) {
  check_mermod_convergence_ipcr(x)
}

#' @noRd
#' @export
check_model_convergence_ipcr.glmerMod <- function(x, ...) {
  check_mermod_convergence_ipcr(x)
}

#' @noRd
#' @export
check_model_convergence_ipcr.MxModel <- function(x, ...) {
  status_code <- x$output$status$code
  if (length(status_code) != 1L || is.na(status_code) || status_code != 0) {
    status_label <- if (length(status_code) == 1L && !is.na(status_code)) {
      as.character(status_code)
    } else {
      "unavailable"
    }
    stop(
      "The fitted OpenMx model did not converge successfully (status code ",
      status_label,
      ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Obtain the number of observations used to compute the casewise scores.
nobs_ipcr <- function(x, ...) {
  UseMethod("nobs_ipcr")
}

#' @noRd
#' @export
nobs_ipcr.default <- function(x, ...) {
  stats::nobs(x, ...)
}

#' @noRd
#' @export
nobs_ipcr.lavaan <- function(x, ...) {
  group_n <- tryCatch(
    lavaan::lavInspect(x, what = "nobs"),
    error = function(e) {
      stop(
        "Unable to obtain the observation count from the fitted lavaan model: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  as.integer(sum(unlist(group_n, use.names = FALSE)))
}

#' @noRd
#' @export
nobs_ipcr.lmerMod <- function(x, ...) {
  nlevels(mixed_model_group_ipcr(x))
}

#' @noRd
#' @export
nobs_ipcr.glmerMod <- function(x, ...) {
  nlevels(mixed_model_group_ipcr(x))
}

#' @noRd
#' @export
nobs_ipcr.MxModel <- function(x, ...) {
  as.integer(NROW(openmx_raw_data_ipcr(x)))
}


validate_scores_ipcr <- function(param_estimates, scores, n) {
  if (!is.numeric(param_estimates) || length(param_estimates) == 0L) {
    stop(
      "The fitted model must provide a non-empty numeric parameter vector.",
      call. = FALSE
    )
  }
  if (any(!is.finite(param_estimates))) {
    stop("The parameter vector contains non-finite values.", call. = FALSE)
  }

  param_names <- names(param_estimates)
  if (is.null(param_names) || anyNA(param_names) || any(!nzchar(param_names))) {
    stop(
      "All fitted model parameters must have non-empty names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(param_names)) {
    stop("Fitted model parameter names must be unique.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n %% 1 != 0) {
    stop("The fitted model must report one positive observation count.", call. = FALSE)
  }

  q <- length(param_estimates)
  scores <- as.matrix(scores)

  if (!is.numeric(scores)) {
    stop("The score matrix must be numeric.", call. = FALSE)
  }
  if (NROW(scores) != n) {
    stop(
      "The score matrix must have ",
      n,
      " rows, one for each IPC unit; it has ",
      NROW(scores),
      ".",
      call. = FALSE
    )
  }
  if (NCOL(scores) != q) {
    stop(
      "The score matrix must have ",
      q,
      " columns, one for each fitted model parameter; it has ",
      NCOL(scores),
      ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(scores))) {
    stop("The score matrix contains non-finite values.", call. = FALSE)
  }

  score_names <- colnames(scores)
  has_complete_names <- !is.null(score_names) &&
    length(score_names) == q &&
    !anyNA(score_names) &&
    all(nzchar(score_names))
  if (has_complete_names && !identical(score_names, param_names)) {
    stop(
      "Score column names do not match the fitted model parameter names or order.",
      call. = FALSE
    )
  }

  colnames(scores) <- param_names
  scores
}


# Validate and normalize the parameter, score, and bread objects used to
# construct IPCs. Complete dimension names are checked when a backend supplies
# them; missing names are filled from the fitted parameter vector.
validate_ipcr_components <- function(param_estimates, scores, bread_matrix, n) {
  scores <- validate_scores_ipcr(
    param_estimates = param_estimates,
    scores = scores,
    n = n
  )
  param_names <- names(param_estimates)
  q <- length(param_estimates)
  bread_matrix <- as.matrix(bread_matrix)

  if (!is.numeric(bread_matrix)) {
    stop("The bread matrix must be numeric.", call. = FALSE)
  }
  if (!identical(dim(bread_matrix), c(q, q))) {
    stop(
      "The bread matrix must be a ",
      q,
      " by ",
      q,
      " matrix; its dimensions are ",
      paste(dim(bread_matrix), collapse = " by "),
      ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(bread_matrix))) {
    stop("The bread matrix contains non-finite values.", call. = FALSE)
  }

  check_names <- function(component_names, component) {
    has_complete_names <- !is.null(component_names) &&
      length(component_names) == q &&
      !anyNA(component_names) &&
      all(nzchar(component_names))

    if (has_complete_names && !identical(component_names, param_names)) {
      stop(
        component,
        " names do not match the fitted model parameter names or order.",
        call. = FALSE
      )
    }
  }

  check_names(colnames(bread_matrix), "Bread column")
  check_names(rownames(bread_matrix), "Bread row")

  if (qr(bread_matrix)$rank < q) {
    stop("The bread matrix is singular.", call. = FALSE)
  }

  dimnames(bread_matrix) <- list(param_names, param_names)

  list(
    param_estimates = param_estimates,
    scores = scores,
    bread_matrix = bread_matrix,
    n = as.integer(n),
    q = q
  )
}
