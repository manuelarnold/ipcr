#' Iterated Individual Parameter Contribution Regression
#'
#' Performs iterated individual parameter contribution regression. The current
#' implementation supports continuous ML models fitted with lavaan and
#' raw-data OpenMx models. OpenMx RAM models can use fast analytical scores;
#' other OpenMx models use a substantially slower numerical fallback.
#'
#' @param x A fitted model object.
#' @param predictors A vector, matrix, or `data.frame` containing one or more
#'   predictors of parameter heterogeneity. Rows must correspond, in order, to
#'   the raw-data rows used by `x`.
#' @param analytic A single logical value. For `MxRAMModel` objects, `TRUE`
#'   uses the fast analytical score calculation and `FALSE` uses repeated
#'   numerical OpenMx score and Jacobian calculations. The numerical route can
#'   be extremely slow. Unlike standard `ipcr()`, analytical iteration rejects
#'   RAM models containing unsupported analytical features rather than silently
#'   switching algorithms. General `MxModel` objects always use numerical row
#'   gradients, regardless of this argument. The lavaan method uses lavaan's
#'   scores for the initial IPCR and analytical RAM calculations during
#'   iteration.
#' @param conv A positive numeric convergence tolerance. An attempt converges
#'   when no second-stage coefficient changes by more than `conv` times its
#'   standard-IPCR standard error after accounting for the learning rate. Thus,
#'   the default `conv = 0.01` limits changes to 1 percent of their initial
#'   standard errors.
#' @param learning_rate A non-empty numeric vector of learning rates in `(0, 1]`.
#'   The default is `c(1, 0.5, 0.1)`. Rates are tried in the supplied order.
#'   Each rate is held fixed during one attempt, and a new attempt restarts from
#'   the pooled fitted model and its score calculation. The function returns as
#'   soon as an attempt converges.
#' @param max_it A positive whole number giving the maximum iterations for each
#'   learning-rate attempt.
#' @param iteration_info A single logical value. If `TRUE`, retain the complete
#'   iteration history in a successful result or in the specialized error
#'   condition signaled when no learning-rate attempt converges.
#' @param variance_parameterization Character. The default, `"original"`,
#'   iterates all model parameters on their fitted scale. For an analytical
#'   OpenMx `MxRAMModel`, `"log"` identifies free parameters that occur
#'   exclusively on the diagonal of the RAM S matrix and iterates them on the
#'   log-variance scale. This option is not available for lavaan, general
#'   `MxModel`, or numerical `MxRAMModel` iteration.
#' @details
#' For observation \eqn{i}, iteration \eqn{t}, and learning rate \eqn{alpha},
#' the recalculated contribution is
#' \deqn{IPC_i^{(t)} = theta_i^{(t)} + alpha I(theta_i^{(t)})^{-1}
#' S_i(theta_i^{(t)}).}
#' The parameter vector \eqn{theta_i^{(t)}} is predicted by the current IPC
#' regression. Smaller learning rates damp the Fisher-scoring correction. If an
#' attempt fails numerically or reaches `max_it`, the next supplied rate is
#' tried. The first converged attempt is returned. If no attempt converges, all
#' supplied rates are exhausted and the function signals an error rather than
#' returning coefficients or inferential results from a nonconverged iteration.
#'
#' For coefficient \eqn{j}, the convergence statistic is
#' \deqn{C_{SE}^{(t)} = \max_j
#'   \frac{|\beta_j^{(t)} - \beta_j^{(t-1)}|}
#'   {\alpha SE(\beta_j^{standard})}.}
#' Dividing by the learning rate \eqn{\alpha} prevents a strongly damped update
#' from being declared converged merely because its step was small. All
#' standard-IPCR coefficient standard errors, including intercept standard
#' errors, must be positive and finite. The diagnostic
#' \deqn{C_{l}^{(t)} = \frac{|l^{(t)} - l^{(t-1)}|}{n}}
#' records the absolute change in rowwise log-likelihood per analyzed
#' observation. It is not part of the stopping rule. When \eqn{C_{SE}} meets
#' the requested convergence tolerance but the learning-rate-adjusted value
#' \eqn{C_l / \alpha} exceeds 0.001, the result is returned with a warning that
#' the likelihood has not stabilized. This is a numerical diagnostic rather
#' than a statistically calibrated test. The `print()` and `summary()` methods
#' repeat this warning so that it remains visible when the result is inspected.
#'
#' Predictor rows with missing values are excluded from the second-stage
#' regression and from every iterative update. Their rows remain in the
#' returned `IPCs` and `scores` matrices, where they retain the values computed
#' by the initial standard IPCR analysis. The `complete_rows` component
#' identifies the rows that were actually iterated.
#'
#' With `variance_parameterization = "log"`, fixed S-diagonal elements remain
#' unchanged and shared labels across multiple S-diagonal elements remain one
#' equality-constrained parameter. A free label that also occurs off the S
#' diagonal or in another RAM matrix is rejected because it cannot
#' simultaneously represent a log variance and an untransformed parameter.
#' The primary IPCs, scores, coefficients, tests, and convergence criterion are
#' reported in the mixed coordinate system: log scale for identified variance
#' parameters and original scale for all other parameters. The returned
#' `natural_scale_IPCs` matrix exponentiates the log-variance IPC columns.
#' The supplementary `output$natural_scale_effects` table reports variance
#' ratios and changes from a one-unit increase in each second-stage design
#' column, with the other design columns fixed at zero. Positive diagonal
#' elements alone do not guarantee that a covariance matrix is positive
#' definite, so log parameterization cannot prevent every numerical failure.
#'
#' The lavaan method currently supports single-group, continuous models fitted
#' by normal-theory ML. Mean structures, simple equality constraints, and FIML
#' data are supported. Multigroup, categorical, non-ML, conditional-x, and
#' sampling-weighted lavaan models are rejected. The OpenMx analytical
#' restrictions are described in [ipcr()]; setting `analytic = FALSE` provides
#' a much slower numerical fallback for `MxRAMModel` objects.
#'
#' For a non-RAM `MxModel`, the numerical method estimates the information at
#' each distinct predicted parameter vector from the outer product of OpenMx's
#' numerical row gradients. It is restricted to converged models with
#' top-level, unweighted raw data and a meaningful row-additive likelihood.
#' Explicit `mxConstraint()` objects and penalties are not supported. The
#' method evaluates the likelihood for each corresponding subset of raw-data
#' rows when validating iterations and recording diagnostics. Models fitted
#' only to summary statistics, models with unusable row gradients, and models
#' whose free parameters cannot be set at the predicted values are rejected. A
#' continuous predictor may create a distinct parameter vector for nearly every row and
#' can consequently make this fallback prohibitively slow.
#'
#' @return On convergence, an object inheriting from classes `"ipcr_it"` and
#'   `"ipcr"`. Its `convergence` component records the selected learning rate
#'   and iteration, the final \eqn{C_{SE}}, \eqn{C_l}, and rate-adjusted
#'   \eqn{C_l / \alpha} diagnostics, a `likelihood_stable` flag, and a summary
#'   of every rate attempted before convergence. If no supplied rate
#'   converges, the function signals an error condition inheriting from
#'   `"ipcr_it_nonconvergence"`. The condition contains an `attempts` data
#'   frame, and also an `iteration_history` data frame when
#'   `iteration_info = TRUE`; it contains no IPC regression coefficients,
#'   tests, or ordinary `"ipcr_it"` result. With log-variance parameterization,
#'   a successful object additionally contains
#'   `natural_scale_IPCs`, `parameterization`, and the supplementary table
#'   `output$natural_scale_effects`. The latter reports the exponentiated
#'   coefficient as a variance ratio and the natural-scale change from setting
#'   the corresponding second-stage design column from zero to one while all
#'   other non-intercept design columns are held at zero.
#'
#' @references
#' Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2020).
#' Identifying heterogeneity in dynamic panel models with individual parameter
#' contribution regression. *Structural Equation Modeling, 27*, 613--628.
#'
#' @examples
#' # A compact iterated lavaan analysis with a binary IPC predictor
#' lavaan_data <- lavaan::HolzingerSwineford1939
#' lavaan_fit <- lavaan::cfa(
#'   "visual =~ x1 + x2 + x3",
#'   data = lavaan_data
#' )
#' lavaan_ipcr_it <- ipcr_it(
#'   lavaan_fit,
#'   predictors = data.frame(sex = lavaan_data$sex - 1)
#' )
#' summary(lavaan_ipcr_it)
#'
#' @export
ipcr_it <- function(
    x, predictors, analytic = TRUE, conv = 0.01,
    learning_rate = c(1, 0.5, 0.1),
    max_it = 350L, iteration_info = FALSE,
    variance_parameterization = "original") {
  UseMethod("ipcr_it")
}

validate_original_variance_parameterization_ipcr <- function(
    variance_parameterization, backend) {
  if (!is.character(variance_parameterization) ||
      length(variance_parameterization) != 1L ||
      is.na(variance_parameterization) ||
      !identical(variance_parameterization, "original")) {
    stop(
      "Log-variance parameterization is available only for analytical ",
      "OpenMx MxRAMModel iteration; ", backend,
      " supports only variance_parameterization = \"original\".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
#' @export
ipcr_it.default <- function(
    x, predictors, analytic = TRUE, conv = 0.01,
    learning_rate = c(1, 0.5, 0.1), max_it = 350L,
    iteration_info = FALSE, variance_parameterization = "original") {
  validate_original_variance_parameterization_ipcr(
    variance_parameterization,
    "this model backend"
  )
  stop(
    "Iterated IPCR is currently implemented only for lavaan and OpenMx MxModel objects.",
    call. = FALSE
  )
}
