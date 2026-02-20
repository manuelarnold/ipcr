#' Individual Parameter Contribution Regression
#'
#' Performs individual parameter contribution regression (IPCR), a method for
#' modeling parameter heterogeneity by regressing parameter estimates on
#' predictors. \code{ipcr} is designed for structural equation models fitted
#' with \pkg{lavaan} and \pkg{OpenMx}, but it also supports
#' \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, and
#' \code{\link[lme4]{lmer}} models.
#'
#' @param fit A fitted model object. Supported models include those estimated
#' using \pkg{lavaan}, \pkg{OpenMx}, \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, and \code{\link[lme4]{lmer}}.
#' @param predictors A vector, matrix, or \code{data.frame} containing one or
#' more predictors used to predict variation in model parameters. Interaction
#' and polynomial terms can be included as new variables, which may require
#' centering.Ensure categorical variables are properly coded as factors or dummy
#' variables.
#' @param analytic Logical. If \code{FALSE} (default), functions of
#' \pkg{lavaan}, \pkg{OpenMx}, or \pkg{sandwich} will be used to compute scores.
#' If \code{TRUE}, custom functions will be used. This is only relevant for
#' models fitted with \pkg{OpenMx} where the computation of the scores can take
#' time. Supports \code{MxRAMModel} without algebras.

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
#' IPCs are known to be slightly biased. This bias can be corrected using
#' iterated IPC regression, which iteratively recalculates the IPCs. The
#' function \code{ipcr_it} implements iterated IPC regression.
#'
#' To identify the most important predictors of heterogeneity, regularization
#' can be applied by passing the output of \code{ipcr} or \code{ipcr_it}
#' to the \code{ipcr_reg} function.
#'
#' @return
#' An object of class \code{"ipcr"}, which is a list containing the following
#' elements:
#'
#' \tabular{ll}{
#' \code{info} \tab A list with metadata about the \code{ipcr} function call. \cr
#' \code{IPCs} \tab A \code{data.frame} containing individual parameter
#' contributions. \cr
#' \code{mlm} \tab An object of class \code{lm} (for a single parameter) or
#' \code{mlm} (for multiple parameters), representing the regression models
#' fitted for each parameter. \cr
#' \code{mancova} \tab An object of class \code{manova} (if the model contains
#' multiple parameters) or \code{NULL} (if the model contains only a single
#' prameter). \cr
#' \code{output} \tab A list of formatted output tables, which can be examined
#' using \code{print} and \code{summary}. \cr
#' }
#'
#' The function \code{summary} provides an overview of the IPC regression
#' models. Use \code{summary.ipcr} to extract specific information from the
#' model results. The \code{print} function displays the arguments specified in
#' the \code{ipcr} call. The \code{plot} function visualizes the correlations
#' between IPCs and predictors as a heatmap.
#'
#' Several generic functions can be applied to \code{ipcr} objects to extract
#' information from the mlm object, including: \code{coef}, \code{confint},
#' \code{effects}, \code{fitted}, \code{nobs}, \code{predict}, \code{residuals},
#' \code{sigma}, and \code{vcov}.
#'
#' Heteroskedasticity-robust IPC regression can be performed using \code{coeftest} and \code{coefci}
#' from the \pkg{lmtest} package.
#'
#' @examples
#' # Load the Holzinger and Swineford (1939) dataset from the lavaan package
#' HS_data <- lavaan::HolzingerSwineford1939
#'
#' # Remove observation 301 because one of its predictor values is missing
#' HS_data <- HS_data[stats::complete.cases(HS_data), ]
#'
#' # Define a confirmatory factor analysis (CFA) model using lavaan syntax
#' # This model specifies three latent factors: visual, textual, and speed
#' m <- 'visual  =~ x1 + x2 + x3
#'       textual =~ x4 + x5 + x6
#'       speed   =~ x7 + x8 + x9'
#'
#' # Fit the CFA model
#' fit <- lavaan::cfa(model = m, data = HS_data)
#'
#' # Preprocess the predictor variables
#'
#' # Combine the year part (ageyr) and the month part (agemo) of age into a single
#' # variable that measures age in years
#' HS_data$age <- HS_data$ageyr + HS_data$agemo / 12
#'
#' # Convert sex and grade into dummy variables (coded as 0 and 1)
#' HS_data$sex <- HS_data$sex - 1
#' HS_data$grade <- HS_data$grade - 7
#'
#' # Create a new data.frame with the predictor variables
#' predictors <- HS_data[, c("sex", "age", "school", "grade")]
#'
#' # Perform Individual Parameter Contribution Regression (IPCR)
#' res <- ipcr(fit = fit, predictors = predictors)
#'
#' # Plot a heatmap showing correlations between parameters and predictors
#' plot(res)
#'
#' # Display a summary of the IPC regression results
#' summary(res)

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

ipcr <- function(fit, predictors, analytic = FALSE) {

  # Checks ----

  ## Check arguments
  check_arguments_ipcr(fit = fit, predictors = predictors,
                       analytic = analytic)


  # Preprocess predictors ----
  ## Transform predictors into data.frame
  predictors <- as.data.frame(predictors)

  ## Name unnamed predictors
  if (is.null(names(predictors)) | any(is.na(names(predictors)))) {
    warning("Some predictor are not named. Renaming all predictors using the
            order of the data.frame.")
    pred_names <- paste0("predictor", seq_len(NCOL(predictors)))
    colnames(predictors) <- pred_names
  } else {
    pred_names <- names(predictors)
  }


  # Storing object for output ----

  ## Information from the model
  param_estimates <- coef_ipcr(fit)
  param_names <- names(param_estimates)
  n <- nobs(fit)
  q <- length(param_estimates)

  ## ipcr object
  IPCR <- list("info" = list(ipcr_type = "standard",
                             name = deparse(substitute(fit)),
                             class = class(fit),
                             parameters = param_names,
                             predictors = pred_names,
                             analytic = analytic))


  # Individual parameter contribution regression --------

  ## Compute scores
  scores <- estfun_ipcr(fit, analytic = analytic)
  bread_matrix <- bread_ipcr(fit)
  IPCs <- matrix(param_estimates, nrow = n, ncol = q, byrow = TRUE) +
    scores %*% t(bread_matrix)
  colnames(IPCs) <- param_names
  IPCR$IPCs <- IPCs

  ## Check if IPCs and predictors have the same number of rows
  if (NROW(predictors) != NROW(IPCs)) {
    stop("The number of rows in the IPCs and the predictor matrix do not
           match. Ensure that each row of the predictor data corresponds to the
           same row as the data used to fit the model. If necessary, remove rows
           with missing values.")
  }

  ## Regress IPCs on predictors
  mlm <- lm(IPCs ~ ., data = predictors)
  IPCR$mlm <- mlm

  ## MANCOVA
  if (q > 1) {
    IPCR$mancova <- manova(mlm)
  } else {
    IPCR$mancova <- NULL
  }


  # Prepare output --------

  ## data.frame with coefficients of the IPCR equations
  ### Initialize an empty data.frame for the regression coefficients
  ipcr_table <- data.frame()

  ### Initialize empty data.frame for the F tests
  F_tests_table <- data.frame()

  # Loop through models and extract coefficients
  for (i in seq_len(NCOL(mlm$coefficients))) {
    mlm_summary <- summary(mlm)[[i]]

    # Create a data frame for this model's coefficients
    temp_df_coef <- data.frame(
      Parameter = colnames(mlm$coefficients)[i],
      Predictor = rownames(mlm_summary$coefficients),
      Estimate = mlm_summary$coefficients[, "Estimate"],
      Std_Error = mlm_summary$coefficients[, "Std. Error"],
      t_Value = mlm_summary$coefficients[, "t value"],
      p_Value = mlm_summary$coefficients[, "Pr(>|t|)"],
      row.names = NULL
    )

    # Create a data.frame for this model's F-test
    F_Value <- mlm_summary$fstatistic[1]
    DF1 <- mlm_summary$fstatistic[2]
    DF2 <- mlm_summary$fstatistic[3]
    temp_df_F_test <- data.frame(
      Parameter = colnames(mlm$coefficients)[i],
      Rsquared = mlm_summary$r.squared,
      F_Value = F_Value,
      DF1 = DF1,
      DF2 = DF2,
      p_Value = pf(q = F_Value, df1 = DF1, df2 = DF2, lower.tail = FALSE),
      row.names = NULL
    )

    # Combine with the main tables
    ipcr_table <- rbind(ipcr_table, temp_df_coef)
    F_tests_table <- rbind(F_tests_table, temp_df_F_test)
  }

  ## Table with IPCR coefficients
  IPCR$output$ipcr <- ipcr_table

  ## MANCOVA table
  IPCR$output$mancova <- summary(IPCR$mancova)$stats

  ## F-test table
  IPCR$output$F_tests <- F_tests_table


  # Return IPCR object ----

  ## Assign class
  class(IPCR) <- "ipcr"

  ## Return IPCR
  IPCR

}
