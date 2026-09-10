# Manual comparison of standard and iterated IPCR for an OpenMx path model
#
# Two structural parameters vary across observations:
#
#   b21 = 0.35 + 0.25 * group
#   b32 = 0.45 + 0.20 * continuous
#
# `group` is a balanced 0/1 dummy, `continuous` is standard normal, and `noise`
# is an unrelated standard-normal predictor. The pooled OpenMx model does not
# include these predictors; standard and iterated IPCR are used afterward to
# recover their associations with the model parameters.
#
# Run this script interactively from the ipcr package root. To exercise all
# computations without opening plot devices, set IPCR_MANUAL_PLOTS=false.


# Load the development version of ipcr ------------------------------------

if (!requireNamespace("OpenMx", quietly = TRUE)) {
  stop("This manual check requires the OpenMx package.", call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop(
    "This manual check requires pkgload. Run it from the ipcr package root.",
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


# Generate heterogeneous model data --------------------------------------

set.seed(20260825)

n <- 300L
group <- sample(rep(0:1, each = n / 2L))
continuous <- stats::rnorm(n)
noise <- stats::rnorm(n)

b21_individual <- 0.35 + 0.25 * group
b32_individual <- 0.45 + 0.20 * continuous

x1 <- stats::rnorm(n)
y2 <- b21_individual * x1 + stats::rnorm(n, sd = sqrt(0.70))
y3 <- b32_individual * y2 + stats::rnorm(n, sd = sqrt(0.60))

model_data <- data.frame(x1 = x1, y2 = y2, y3 = y3)
predictors <- data.frame(
  group = group,
  continuous = continuous,
  noise = noise
)

generating_effects <- data.frame(
  Parameter = c("b21", "b32"),
  Predictor = c("group", "continuous"),
  Intercept = c(0.35, 0.45),
  Effect = c(0.25, 0.20)
)

cat("\nGenerating parameter relationships:\n")
print(generating_effects, row.names = FALSE)


# Fit the pooled OpenMx RAM model -----------------------------------------

pooled_model <- OpenMx::mxModel(
  "SimplePathModel",
  type = "RAM",
  manifestVars = c("x1", "y2", "y3"),
  OpenMx::mxPath(
    from = "x1",
    to = "y2",
    arrows = 1,
    free = TRUE,
    values = 0.35,
    labels = "b21"
  ),
  OpenMx::mxPath(
    from = "y2",
    to = "y3",
    arrows = 1,
    free = TRUE,
    values = 0.45,
    labels = "b32"
  ),
  OpenMx::mxPath(
    from = "x1",
    arrows = 2,
    free = TRUE,
    values = 1,
    labels = "var_x1",
    lbound = 0.0001
  ),
  OpenMx::mxPath(
    from = "y2",
    arrows = 2,
    free = TRUE,
    values = 0.70,
    labels = "res_y2",
    lbound = 0.0001
  ),
  OpenMx::mxPath(
    from = "y3",
    arrows = 2,
    free = TRUE,
    values = 0.60,
    labels = "res_y3",
    lbound = 0.0001
  ),
  OpenMx::mxPath(
    from = "one",
    to = c("x1", "y2", "y3"),
    arrows = 1,
    free = FALSE,
    values = 0
  ),
  OpenMx::mxData(observed = model_data, type = "raw")
)

fit <- NULL
invisible(utils::capture.output(
  fit <- OpenMx::mxTryHard(
    pooled_model,
    extraTries = 5,
    bestInitsOutput = FALSE,
    silent = TRUE,
    verbose = 0
  ),
  type = "output"
))

cat("\nPooled OpenMx estimates:\n")
print(OpenMx::omxGetParameters(fit))


# Standard and iterated IPCR ---------------------------------------------

standard_ipcr <- ipcr(
  x = fit,
  predictors = predictors,
  analytic = TRUE
)

iterated_ipcr_result <- ipcr_it(
  x = fit,
  predictors = predictors,
  analytic = TRUE,
  iteration_info = TRUE
)

cat("\n\nSTANDARD IPCR SUMMARY\n")
cat("=====================\n\n")
summary(standard_ipcr, what = "everything")

cat("\n\nITERATED IPCR SUMMARY\n")
cat("=====================\n\n")
summary(iterated_ipcr_result, what = "everything")

cat("\n\nITERATED IPCR CONVERGENCE\n")
cat("==========================\n\n")
cat(iterated_ipcr_result$convergence$status, "\n\n")
print(iterated_ipcr_result$convergence$attempts, row.names = FALSE)

cat("\nSelected effects for comparison:\n")
effect_rows <- function(result) {
  result$output$ipcr[
    (result$output$ipcr$Parameter == "b21" &
       result$output$ipcr$Predictor %in% c("group", "continuous", "noise")) |
      (result$output$ipcr$Parameter == "b32" &
         result$output$ipcr$Predictor %in% c("group", "continuous", "noise")),
    c("Parameter", "Predictor", "Estimate", "Std_Error", "p_Value"),
    drop = FALSE
  ]
}

cat("\nStandard IPCR:\n")
print(effect_rows(standard_ipcr), row.names = FALSE)
cat("\nIterated IPCR:\n")
print(effect_rows(iterated_ipcr_result), row.names = FALSE)


# Plot the results --------------------------------------------------------

show_plots <- !identical(
  tolower(Sys.getenv("IPCR_MANUAL_PLOTS", unset = "true")),
  "false"
)

standard_heatmap <- plot(standard_ipcr, print_corr = TRUE)
iterated_heatmap <- plot(iterated_ipcr_result, print_corr = TRUE)

standard_difference_plots <- plot_differences(
  standard_ipcr,
  parameter = c("b21", "b32"),
  predictor = c("group", "continuous", "noise")
)
iterated_difference_plots <- plot_differences(
  iterated_ipcr_result,
  parameter = c("b21", "b32"),
  predictor = c("group", "continuous", "noise")
)

if (show_plots) {
  print(standard_heatmap)
  print(iterated_heatmap)

  for (plot_name in names(standard_difference_plots)) {
    message("Standard IPCR: ", plot_name)
    print(standard_difference_plots[[plot_name]])
  }
  for (plot_name in names(iterated_difference_plots)) {
    message("Iterated IPCR: ", plot_name)
    print(iterated_difference_plots[[plot_name]])
  }
} else {
  message("Plot construction completed; display skipped by IPCR_MANUAL_PLOTS=false.")
}
