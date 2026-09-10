#' Individual Parameter Contribution Regression
#'
#' Performs individual parameter contribution regression (IPCR), a method for
#' modeling parameter heterogeneity by regressing parameter estimates on
#' predictors. The primary backends are \pkg{OpenMx} and \pkg{lavaan}.
#' Additional official backends are
#' \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, and the
#' \code{\link[lme4]{lmer}} and \code{\link[lme4]{glmer}} model classes.
#'
#' @param x A fitted model object. Supported models include those estimated
#' using \pkg{OpenMx}, \pkg{lavaan}, \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, \code{\link[lme4]{lmer}}, and
#' \code{\link[lme4]{glmer}}. Other model classes are attempted through
#' compatible \pkg{sandwich} methods with a warning.
#' @param predictors A numeric, logical, or factor vector, or a matrix or
#'   \code{data.frame} containing one or more such predictors. Character
#'   predictors must first be converted to factors or numeric indicators.
#'   Interaction and polynomial terms can be supplied as additional columns.
#'   For mixed models, predictors must contain one row per level of the single
#'   grouping factor, ordered according to that factor's levels.
#' @param analytic Logical. If \code{FALSE} (default), functions of
#' \pkg{lavaan}, \pkg{OpenMx}, or \pkg{sandwich} will be used to compute scores.
#' If \code{TRUE}, fast analytical scores are requested for a continuous
#' raw-data \code{MxRAMModel}. See Details for the supported analytical scope.
#' In standard IPCR, unsupported OpenMx models use numerical OpenMx row
#' gradients with a warning.

#' @details
#' Individual parameter contributions (IPCs) provide rough approximations of
#' individual-specific parameter values. The IPCs of individual \eqn{i} is
#' defined as:
#' \deqn{IPC_i = \theta + A(\theta)^{-1} S(\theta, y_i),}
#' where \eqn{\theta} represents the estimated model parameters,
#' \eqn{S(\theta, y_i)} is the estimating function (e.g., the first derivative
#' of the log-likelihood), and \eqn{A(\theta)} is the expected negative
#' derivative of the estimating function (often called the Hessian matrix).
#' By regressing IPCs on predictors, the sources of parameter heterogeneity can
#' be examined.
#'
#' For ordinary univariate \code{lm} models, IPCs include regression
#' coefficients and the maximum-likelihood residual variance
#' \eqn{RSS / n}. For \code{glm} models, only the regression coefficients are
#' included; no residual-dispersion parameter is added.
#'
#' Mixed-model IPCs are currently defined at the cluster level. Official
#' \code{lmer} support requires maximum-likelihood estimation
#' (\code{REML = FALSE}) and exactly one grouping factor. Official
#' \code{glmer} support is limited to binomial and Poisson models with exactly
#' one grouping factor. Prior weights and \code{cbind} responses are not yet
#' supported for mixed models.
#'
#' Lavaan scores and bread matrices are obtained entirely from \pkg{lavaan}.
#' Supported continuous-data estimators are those represented internally by
#' lavaan as ML, GLS, ULS, or WLS. This includes multigroup models and
#' full-information ML with incomplete data. Simple linear equality constraints
#' are supported, whether specified through shared labels or explicit equality
#' syntax. Nonlinear and inequality constraints are rejected. Ordered-indicator
#' models are currently supported only with \code{estimator = "WLS"}; lavaan
#' does not currently return compatible casewise scores for the default WLSMV
#' estimator or other estimators represented internally as DWLS. Estimators such
#' as PML that do not provide expected information are also rejected.
#'
#' OpenMx support requires a converged model containing raw data. RAM models use
#' numerical OpenMx row gradients by default. Non-RAM \code{MxModel} objects
#' also use numerical row gradients. Because IPCR requires one score per IPC
#' unit, an OpenMx model must provide a meaningful row-wise likelihood; models
#' fitted only to summary statistics may not provide usable row gradients.
#'
#' Fast analytical OpenMx scores are available for continuous raw-data RAM
#' models whose free parameters occur in the A, S, or M matrices. They support
#' missing manifest observations, definition variables in A, S, and M, and
#' equality restrictions represented by shared parameter labels. Definition
#' variables themselves must be numeric and complete. Models containing
#' explicit \code{mxConstraint()} objects, algebras, free F-matrix parameters,
#' submodels, penalties, ordinal outcomes, or row/frequency weights fall back to
#' numerical scores in \code{ipcr()}, with a warning. Analytical
#' \code{ipcr_it()} instead rejects unsupported RAM features; its explicitly
#' requested \code{analytic = FALSE} route provides numerical iteration.
#' Algebras unrelated to the RAM expectation are also treated as unsupported.
#'
#' Other fitted-model classes are not officially supported. The function will
#' warn and then attempt to use their \pkg{sandwich} score and bread methods.
#' Successful extraction still requires coefficients, scores, and bread to
#' describe the same parameters in the same order. Experimental classes
#' include \code{clm}, \code{coxph}, \code{hurdle}, \code{mlogit},
#' \code{mlm}, \code{nls}, \code{polr}, \code{rlm}, \code{survreg}, and
#' \code{zeroinfl}.
#'
#' Predictor rows are matched positionally to IPC units. Predictor columns must
#' have unique names; unnamed columns are assigned names beginning with
#' \code{predictor1}. If the first-stage model omits observations, predictors
#' must contain only the retained IPC units, in the order used by the fitted
#' model. Rows containing missing predictor values are omitted from the
#' second-stage IPC regression. The complete set of IPCs and scores is
#' nevertheless retained in the returned object.
#'
#' The second-stage coefficient table and overall F-tests are obtained from the
#' separate linear regression for each IPC response. The reported R-squared
#' value describes the proportion of variation in that IPC explained jointly
#' by the predictors, and the F-test tests all non-intercept coefficients for
#' that IPC simultaneously. When there is more than one IPC response and the
#' multivariate test is estimable, \code{\link[car]{Manova}} supplies
#' order-invariant Type-II multivariate tests using Pillai's trace. Each
#' predictor term is therefore tested conditionally on the remaining terms
#' according to the Type-II principle of marginality. The MANCOVA summary is
#' unavailable when the residual IPC response matrix is rank deficient or the
#' multivariate calculation otherwise fails; the separate IPC regressions and
#' their F-tests remain available.
#' These inferential summaries treat the estimated IPCs as generated responses
#' and do not separately propagate first-step estimation uncertainty.
#'
#' IPCs are known to be slightly biased. Iterated IPC regression is intended to
#' reduce this bias by iteratively recalculating the IPCs. See
#' \code{\link{ipcr_it}} for the currently supported lavaan and OpenMx models.
#'
#' Regularized IPCR applies a cross-validated multivariate group elastic net to
#' an existing standard or iterated IPCR result. See \code{\link{ipcr_reg}}.
#'
#' @return
#' An object of class \code{"ipcr"}, which is a list containing the following
#' elements:
#'
#' \tabular{ll}{
#' \code{info} \tab A list with metadata about the \code{ipcr} function call. \cr
#' \code{IPCs} \tab A numeric matrix containing individual parameter
#' contributions. \cr
#' \code{scores} \tab A matrix containing the scores (row-wise gradients of
#' the estimation function). \cr
#' \code{predictors} \tab The full normalized predictor \code{data.frame},
#' including rows with missing predictor values. \cr
#' \code{complete_rows} \tab A logical vector indicating which predictor rows
#' were used in the second-stage regression. \cr
#' \code{mlm} \tab An object of class \code{lm} (for a single parameter) or
#' \code{mlm} (for multiple parameters), representing the regression models
#' fitted for each parameter. \cr
#' \code{mancova} \tab A Type-II multivariate analysis object produced by
#' \code{\link[car]{Manova}} when fitting it succeeds and multiple parameters
#' are present; otherwise \code{NULL}. Its printable summary can be unavailable
#' for a rank-deficient residual IPC response matrix. \cr
#' \code{output} \tab A list of formatted output tables, which can be examined
#' using \code{print} and \code{summary}. \cr
#' }
#'
#' The \code{print} method displays concise model metadata, including the model
#' name and class and the parameter and predictor names. The \code{summary}
#' method prints the marginal F-tests and available MANCOVA, the IPC regression
#' coefficient table, or both, according to its \code{what} argument. The
#' \code{plot} method visualizes the correlations between IPCs and predictors
#' as a heatmap.
#'
#' Several generic functions can be applied to \code{ipcr} objects to extract
#' information from the mlm object, including: \code{coef}, \code{confint},
#' \code{effects}, \code{fitted}, \code{nobs}, \code{predict}, \code{residuals},
#' \code{sigma}, and \code{vcov}.
#'
#' @examples
#' # OpenMx: a raw-data RAM model with analytical scores
#' set.seed(1)
#' n <- 60L
#' openmx_data <- data.frame(
#'   y = stats::rnorm(n),
#'   group = rep(c(0, 1), each = n / 2)
#' )
#' openmx_model <- OpenMx::mxModel(
#'   "one_variable",
#'   type = "RAM",
#'   manifestVars = "y",
#'   OpenMx::mxPath(
#'     from = "y", arrows = 2, free = TRUE, values = 1,
#'     labels = "variance", lbound = 0.0001
#'   ),
#'   OpenMx::mxPath(
#'     from = "one", to = "y", arrows = 1, free = TRUE,
#'     values = 0, labels = "mean"
#'   ),
#'   OpenMx::mxData(openmx_data["y"], type = "raw")
#' )
#' openmx_fit <- OpenMx::mxTryHard(openmx_model, extraTries = 3)
#' openmx_ipcr <- ipcr(
#'   openmx_fit,
#'   predictors = openmx_data["group"],
#'   analytic = TRUE
#' )
#'
#' # lavaan: a small confirmatory factor model
#' lavaan_data <- lavaan::HolzingerSwineford1939
#' lavaan_fit <- lavaan::cfa(
#'   "visual =~ x1 + x2 + x3",
#'   data = lavaan_data
#' )
#' lavaan_ipcr <- ipcr(
#'   lavaan_fit,
#'   predictors = data.frame(sex = lavaan_data$sex - 1)
#' )
#'
#' # lm: slope heterogeneity by a group predictor
#' set.seed(2)
#' n <- 100L
#' group <- rep(c(0, 1), each = n / 2)
#' x <- stats::rnorm(n)
#' y <- 1 + 0.5 * x * group + stats::rnorm(n)
#' lm_fit <- stats::lm(y ~ x)
#' lm_ipcr <- ipcr(lm_fit, predictors = data.frame(group = group))

#' @references
#' Arnold, M., Oberski, D. L., Brandmaier, A. M., & Voelkle, M. C. (2020).
#' Identifying heterogeneity in dynamic panel models with individual parameter
#' contribution regression. \emph{Structural Equation Modeling, 27}, 613-628.
#' doi: \href{https://doi.org/10.1080/10705511.2019.1667240}{10.1080/10705511.2019.1667240}
#'
#' Arnold, M., Brandmaier, A. M., & Voelkle, M. C. (2021). Predicting
#' differences in model parameters with individual parameter contribution
#' regression using the R package ipcr. \emph{Psych, 3}, 360-385.
#' doi: \href{https://doi.org/10.3390/psych3030027}{10.3390/psych3030027}
#'
#' @seealso \code{\link{ipcr_it}}, \code{\link{ipcr_reg}}
#' @export

ipcr <- function(x, predictors, analytic = FALSE) {

  # Checks ----

  ## Check arguments
  analytic <- check_arguments_ipcr(
    x = x,
    predictors = predictors,
    analytic = analytic
  )
  check_model_convergence_ipcr(x)


  # Storing object for output ----

  ## Information from the model
  n <- nobs_ipcr(x)
  predictor_components <- prepare_predictors_ipcr(predictors, n = n)
  predictors <- predictor_components$predictors
  analysis_predictors <- predictor_components$analysis_predictors
  complete_rows <- predictor_components$complete_rows
  pred_names <- names(predictors)
  param_estimates <- coef_ipcr(x)
  param_names <- names(param_estimates)
  q <- length(param_estimates)

  ## ipcr object
  IPCR <- list("info" = list(ipcr_type = "standard",
                             name = deparse(substitute(x)),
                             class = class(x),
                             parameters = param_names,
                             predictors = pred_names,
                             analytic = analytic))


  # Individual parameter contribution regression --------

  ## Compute scores
  scores <- estfun_ipcr(x, analytic = analytic)
  bread_matrix <- bread_ipcr(x)
  ipc_components <- validate_ipcr_components(
    param_estimates = param_estimates,
    scores = scores,
    bread_matrix = bread_matrix,
    n = n
  )
  scores <- ipc_components$scores
  bread_matrix <- ipc_components$bread_matrix
  IPCs <- matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
    scores %*% t(bread_matrix)
  colnames(IPCs) <- param_names
  IPCR$IPCs <- IPCs
  IPCR$scores <- scores
  IPCR$predictors <- predictors
  IPCR$complete_rows <- complete_rows

  ## Regress IPCs on predictors
  analysis_IPCs <- IPCs[complete_rows, , drop = FALSE]
  mlm <- fit_second_stage_ipcr(
    IPCs = analysis_IPCs,
    predictors = analysis_predictors
  )
  IPCR$mlm <- mlm

  # Prepare output --------
  second_stage <- second_stage_output_ipcr(
    mlm = mlm,
    parameter_names = param_names
  )
  IPCR$mancova <- second_stage$mancova
  IPCR$output <- second_stage$output


  # Return IPCR object ----

  ## Assign class
  class(IPCR) <- "ipcr"

  ## Return IPCR
  IPCR

}
