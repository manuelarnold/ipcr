# Manual comparison of iterated IPCR backends
#
# The same univariate normal location-scale model is fitted as:
#
#   1. an MxRAMModel with analytical scores,
#   2. the same MxRAMModel with numerical scores,
#   3. a non-RAM MxModel with numerical scores, and
#   4. a lavaan model.
#
# Only a balanced binary IPC predictor is used deliberately. A continuous
# predictor can produce a unique predicted parameter vector for every row and
# make the two numerical OpenMx routes take a very long time. This script is a
# manual backend comparison and is not part of the automatic test suite.


# Load the development version of ipcr ------------------------------------

required_packages <- c("OpenMx", "lavaan", "pkgload")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_packages) > 0L) {
  stop(
    "This manual check requires: ",
    paste(missing_packages, collapse = ", "),
    ".",
    call. = FALSE
  )
}
if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the ipcr package root.", call. = FALSE)
}
description <- read.dcf("DESCRIPTION", fields = "Package")
if (!identical(unname(description[1, 1]), "ipcr")) {
  stop("Run this script from the ipcr package root.", call. = FALSE)
}

pkgload::load_all(".", quiet = TRUE)
options(width = 120)


# Generate one heterogeneous normal outcome ------------------------------

set.seed(20260826)
n <- 24L
group <- rep(c(0, 1), each = n / 2L)
data <- data.frame(
  y = 0.35 + 0.55 * group + stats::rnorm(
    n,
    sd = ifelse(group == 0, 0.8, 1.05)
  )
)
predictors <- data.frame(group = group)


# Fit equivalent OpenMx and lavaan models --------------------------------

fit_openmx_manual <- function(model) {
  fit <- NULL
  invisible(utils::capture.output(
    fit <- OpenMx::mxTryHard(
      model,
      extraTries = 3,
      bestInitsOutput = FALSE,
      silent = TRUE,
      verbose = 0
    ),
    type = "output"
  ))
  fit
}

ram_model <- OpenMx::mxModel(
  "RAMLocationScale",
  type = "RAM",
  manifestVars = "y",
  OpenMx::mxPath(
    from = "y",
    arrows = 2,
    free = TRUE,
    values = 1,
    labels = "variance",
    lbound = 0.0001
  ),
  OpenMx::mxPath(
    from = "one",
    to = "y",
    arrows = 1,
    free = TRUE,
    values = 0,
    labels = "mean"
  ),
  OpenMx::mxData(data, type = "raw")
)

nonram_model <- OpenMx::mxModel(
  "NonRAMLocationScale",
  OpenMx::mxMatrix(
    type = "Symm",
    nrow = 1,
    ncol = 1,
    free = TRUE,
    values = 1,
    labels = "variance",
    lbound = 0.0001,
    name = "expected_covariance"
  ),
  OpenMx::mxMatrix(
    type = "Full",
    nrow = 1,
    ncol = 1,
    free = TRUE,
    values = 0,
    labels = "mean",
    name = "expected_means"
  ),
  OpenMx::mxData(data, type = "raw"),
  OpenMx::mxExpectationNormal(
    covariance = "expected_covariance",
    means = "expected_means",
    dimnames = "y"
  ),
  OpenMx::mxFitFunctionML()
)

ram_fit <- fit_openmx_manual(ram_model)
nonram_fit <- fit_openmx_manual(nonram_model)
lavaan_fit <- lavaan::sem(
  model = "
    y ~ 1
    y ~~ y
  ",
  data = data,
  meanstructure = TRUE
)


# Run all four iterated backends -----------------------------------------

run_timed_ipcr_it <- function(x, analytic) {
  result <- NULL
  elapsed <- system.time(
    result <- ipcr_it(
      x = x,
      predictors = predictors,
      analytic = analytic,
      conv = 0.01,
      max_it = 15,
      iteration_info = TRUE
    )
  )
  list(result = result, elapsed = unname(elapsed["elapsed"]))
}

ram_analytic_run <- run_timed_ipcr_it(ram_fit, analytic = TRUE)
ram_numeric_run <- run_timed_ipcr_it(ram_fit, analytic = FALSE)
mxmodel_numeric_run <- run_timed_ipcr_it(nonram_fit, analytic = FALSE)
lavaan_run <- run_timed_ipcr_it(lavaan_fit, analytic = TRUE)

iterated_backend_results <- list(
  RAM_analytic = ram_analytic_run$result,
  RAM_numeric = ram_numeric_run$result,
  MxModel_numeric = mxmodel_numeric_run$result,
  lavaan = lavaan_run$result
)


# Print summaries and compact comparisons --------------------------------

for (backend in names(iterated_backend_results)) {
  cat("\n\n", backend, "\n", sep = "")
  cat(strrep("=", nchar(backend)), "\n", sep = "")
  print(summary(iterated_backend_results[[backend]]))
}

runtime_comparison <- data.frame(
  backend = names(iterated_backend_results),
  elapsed_seconds = c(
    ram_analytic_run$elapsed,
    ram_numeric_run$elapsed,
    mxmodel_numeric_run$elapsed,
    lavaan_run$elapsed
  ),
  converged = vapply(
    iterated_backend_results,
    function(result) isTRUE(result$convergence$converged),
    logical(1)
  ),
  learning_rate = vapply(
    iterated_backend_results,
    function(result) result$convergence$learning_rate,
    numeric(1)
  ),
  iterations = vapply(
    iterated_backend_results,
    function(result) result$info$iterations,
    integer(1)
  ),
  row.names = NULL
)

extract_group_effects <- function(result, variance_name, mean_name) {
  coefficients <- stats::coef(result)
  c(
    variance = unname(coefficients["group", variance_name]),
    mean = unname(coefficients["group", mean_name])
  )
}

effect_comparison <- rbind(
  RAM_analytic = extract_group_effects(
    iterated_backend_results$RAM_analytic,
    "variance",
    "mean"
  ),
  RAM_numeric = extract_group_effects(
    iterated_backend_results$RAM_numeric,
    "variance",
    "mean"
  ),
  MxModel_numeric = extract_group_effects(
    iterated_backend_results$MxModel_numeric,
    "variance",
    "mean"
  ),
  lavaan = extract_group_effects(
    iterated_backend_results$lavaan,
    "y~~y",
    "y~1"
  )
)

cat("\n\nRuntime and convergence comparison\n")
print(runtime_comparison, row.names = FALSE)
cat("\nEstimated group effects\n")
print(effect_comparison)


# Construct and optionally display heatmaps ------------------------------

backend_heatmaps <- lapply(
  iterated_backend_results,
  plot,
  print_corr = TRUE
)
show_plots <- !identical(
  tolower(Sys.getenv("IPCR_MANUAL_PLOTS", unset = "true")),
  "false"
)
if (show_plots) {
  for (backend in names(backend_heatmaps)) {
    message("Iterated IPCR heatmap: ", backend)
    print(backend_heatmaps[[backend]])
  }
} else {
  message("Plot construction completed; display skipped by IPCR_MANUAL_PLOTS=false.")
}
