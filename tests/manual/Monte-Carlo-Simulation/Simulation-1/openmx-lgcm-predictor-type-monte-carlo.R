# Predictor-type Monte Carlo simulation for an OpenMx latent growth curve
# model (LGCM).
#
# The informative predictor is either a balanced 0/1 dummy or U(0, 1). In both
# conditions the residual variance is 0.3 + 0.2 * informative. Two independent
# noise predictors, one balanced dummy and one U(0, 1), are included in every
# IPC regression. Standard and iterated IPCR both use analytical RAM scores.
#
# Run from the package root. For a short smoke run:
# IPCR_MC_REPLICATIONS=2 Rscript tests/manual/Monte-Carlo-Simulation/Simulation-1/openmx-lgcm-predictor-type-monte-carlo.R


# Settings and package loading --------------------------------------------

group_size <- 125L
sample_size <- 2L * group_size
n_replications <- as.integer(Sys.getenv("IPCR_MC_REPLICATIONS", "1000"))
checkpoint_every <- as.integer(Sys.getenv("IPCR_MC_CHECKPOINT_EVERY", "10"))
simulation_seed <- 23484L
simulation_directory <- file.path(
  "tests", "manual", "Monte-Carlo-Simulation", "Simulation-1"
)
results_file <- Sys.getenv(
  "IPCR_MC_RESULTS_FILE",
  file.path(
    simulation_directory, "Results",
    "openmx-lgcm-predictor-type-n250-results.rds"
  )
)

if (is.na(n_replications) || n_replications < 1L ||
    is.na(checkpoint_every) || checkpoint_every < 1L) {
  stop("Replication and checkpoint counts must be positive integers.")
}
required_packages <- c("OpenMx", "pkgload")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
)]
if (length(missing_packages) > 0L) {
  stop("This simulation requires: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(OpenMx))
if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] != "ipcr") {
  stop("Run this script from the ipcr package root.")
}

pkgload::load_all(".", quiet = TRUE)
ipcr_it_defaults <- formals(ipcr_it)
iterated_convergence <- eval(ipcr_it_defaults$conv)
iterated_learning_rate <- eval(ipcr_it_defaults$learning_rate)
iterated_max_iterations <- eval(ipcr_it_defaults$max_it)
source(file.path(simulation_directory, "ipcr-monte-carlo-helpers.R"))


# Population targets and pooled OpenMx model ------------------------------

manifest_variables <- paste0("y", 1:4)
time_scores <- 0:3
target_intercepts <- c(
  residual = 0.3, VarI = 1, CovIS = 0.2,
  VarS = 0.1, MeanI = 1, MeanS = 0.1
)
target_effects <- c(
  residual = 0.2, VarI = 0, CovIS = 0,
  VarS = 0, MeanI = 0, MeanS = 0
)
latent_mean <- c(1, 0.1)
latent_covariance <- matrix(c(1, 0.2, 0.2, 0.1), nrow = 2L)

rmvnorm_one <- function(mean, covariance) {
  mean + as.vector(t(chol(covariance)) %*% stats::rnorm(length(mean)))
}

simulate_lgcm_data <- function(informative) {
  data <- matrix(
    NA_real_, nrow = length(informative), ncol = length(time_scores),
    dimnames = list(NULL, manifest_variables)
  )
  residual_variance <- 0.3 + 0.2 * informative
  for (individual in seq_along(informative)) {
    random_effects <- rmvnorm_one(latent_mean, latent_covariance)
    data[individual, ] <- random_effects[1L] +
      random_effects[2L] * time_scores +
      stats::rnorm(
        length(time_scores),
        sd = sqrt(residual_variance[individual])
      )
  }
  as.data.frame(data)
}

make_sample_model <- function(data) {
  OpenMx::mxModel(
    "PooledLGCM",
    type = "RAM",
    manifestVars = manifest_variables,
    latentVars = c("intercept", "slope"),
    OpenMx::mxPath(
      from = "intercept", to = manifest_variables,
      arrows = 1, free = FALSE, values = 1
    ),
    OpenMx::mxPath(
      from = "slope", to = manifest_variables,
      arrows = 1, free = FALSE, values = time_scores
    ),
    OpenMx::mxPath(
      from = manifest_variables, arrows = 2, free = TRUE,
      values = 0.4, labels = "residual", lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "intercept", arrows = 2, free = TRUE,
      values = 1, labels = "VarI", lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "slope", arrows = 2, free = TRUE,
      values = 0.1, labels = "VarS", lbound = 0.0001
    ),
    OpenMx::mxPath(
      from = "intercept", to = "slope", arrows = 2,
      free = TRUE, values = 0.2, labels = "CovIS"
    ),
    OpenMx::mxPath(
      from = "one", to = c("intercept", "slope"), arrows = 1,
      free = TRUE, values = latent_mean, labels = c("MeanI", "MeanS")
    ),
    OpenMx::mxData(data, type = "raw")
  )
}

fit_openmx <- function(model) {
  fit <- NULL
  invisible(utils::capture.output(
    fit <- suppressMessages(OpenMx::mxTryHard(
      model, extraTries = 5, silent = TRUE, verbose = 0,
      bestInitsOutput = FALSE
    )),
    type = "output"
  ))
  if (inherits(fit, "try-error") || !inherits(fit, "MxModel")) {
    stop("OpenMx did not return a fitted MxModel.")
  }
  status_code <- fit$output$status$code
  if (length(status_code) != 1L || is.na(status_code) || status_code != 0) {
    status_label <- if (length(status_code) == 1L && !is.na(status_code)) {
      as.character(status_code)
    } else {
      "unavailable"
    }
    stop(
      "The OpenMx model did not converge successfully (status code ",
      status_label, ")."
    )
  }
  fit
}


# One replication ---------------------------------------------------------

run_replication <- function(replication) {
  set.seed(simulation_seed + replication - 1L)
  conditions <- list(
    informative_dummy = sample(rep(0:1, each = group_size)),
    informative_continuous = stats::runif(sample_size, 0, 1)
  )
  estimates <- diagnostics <- tests <- failures <- list()
  result_index <- 0L

  for (condition in names(conditions)) {
    informative <- conditions[[condition]]
    predictors <- data.frame(
      informative = informative,
      noise_dummy = sample(rep(0:1, each = group_size)),
      noise_uniform = stats::runif(sample_size, 0, 1)
    )
    model_data <- simulate_lgcm_data(informative)
    fit_attempt <- mc_capture(fit_openmx(make_sample_model(model_data)))
    if (inherits(fit_attempt$value, "error")) {
      failures[[length(failures) + 1L]] <- mc_failure(
        replication, "pooled_model", conditionMessage(fit_attempt$value),
        condition = condition
      )$failures
      next
    }
    fit <- fit_attempt$value
    if (!identical(names(OpenMx::omxGetParameters(fit)), names(target_effects))) {
      stop("Unexpected OpenMx parameter names or order.")
    }

    for (method in c("standard", "iterated")) {
      result_index <- result_index + 1L
      analysis_attempt <- if (method == "standard") {
        mc_capture(suppressMessages(ipcr::ipcr(
          fit, predictors = predictors, analytic = TRUE
        )))
      } else {
        mc_capture(suppressMessages(ipcr::ipcr_it(
          fit,
          predictors = predictors,
          analytic = TRUE,
          iteration_info = FALSE
        )))
      }
      if (inherits(analysis_attempt$value, "error")) {
        failures[[length(failures) + 1L]] <- mc_failure(
          replication, "ipcr_analysis", conditionMessage(analysis_attempt$value),
          condition = condition, method = method,
          error = analysis_attempt$value
        )$failures
        next
      }

      result <- analysis_attempt$value
      estimates[[result_index]] <- mc_extract_ipcr_estimates(
        result, replication, condition, method,
        target_intercepts, target_effects, informative_name = "informative"
      )
      for (predictor in names(predictors)) {
        pillai <- mc_pillai_test(result, predictor)
        tests[[length(tests) + 1L]] <- data.frame(
          replication = replication,
          condition = condition,
          method = method,
          test = paste0("global_", predictor),
          statistic = pillai[["statistic"]],
          df1 = pillai[["df1"]],
          df2 = pillai[["df2"]],
          p_value = pillai[["p_value"]],
          row.names = NULL
        )
      }
      is_iterated <- method == "iterated"
      iteration_diagnostics <- mc_iteration_diagnostics(result, is_iterated)
      diagnostics[[result_index]] <- data.frame(
        replication = replication,
        condition = condition,
        method = method,
        pooled_fit_seconds = fit_attempt$elapsed_seconds,
        elapsed_seconds = analysis_attempt$elapsed_seconds,
        converged = iteration_diagnostics$converged,
        learning_rate = iteration_diagnostics$learning_rate,
        iterations = iteration_diagnostics$iterations,
        total_iterations = iteration_diagnostics$total_iterations,
        learning_rate_attempts = iteration_diagnostics$learning_rate_attempts,
        C_SE = iteration_diagnostics$C_SE,
        C_l = iteration_diagnostics$C_l,
        attempts = I(list(iteration_diagnostics$attempts)),
        pooled_fit_warnings = paste(fit_attempt$warnings, collapse = " | "),
        warnings = paste(analysis_attempt$warnings, collapse = " | "),
        row.names = NULL
      )
    }
  }

  list(
    estimates = mc_combine_rows(estimates),
    diagnostics = mc_combine_rows(diagnostics),
    tests = mc_combine_rows(tests),
    failures = mc_combine_rows(failures)
  )
}


# Simulation run ----------------------------------------------------------

settings <- list(
  design = "OpenMx LGCM predictor-type simulation",
  sample_size = sample_size,
  n_replications = n_replications,
  checkpoint_every = checkpoint_every,
  seed = simulation_seed,
  conditions = c("informative_dummy", "informative_continuous"),
  methods = c("standard", "iterated"),
  informative_continuous_distribution = "Uniform(0, 1)",
  noise_predictors = c("balanced dummy", "Uniform(0, 1)"),
  global_tests = c(
    "Pillai MANCOVA for informative",
    "Pillai MANCOVA for noise_dummy",
    "Pillai MANCOVA for noise_uniform"
  ),
  target_intercepts = target_intercepts,
  target_effects = target_effects,
  analytic = TRUE,
  iterated_convergence = iterated_convergence,
  iterated_learning_rate = iterated_learning_rate,
  iterated_max_iterations = iterated_max_iterations,
  ipcr_version = as.character(utils::packageVersion("ipcr")),
  OpenMx_version = as.character(utils::packageVersion("OpenMx"))
)

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
resume <- mc_resume_checkpoint(results_file, settings, n_replications)
replication_results <- resume$results
progress_counts <- resume$counts
started_at <- Sys.time()
if (resume$completed_replications > 0L) {
  message(
    "Resuming the LGCM predictor-type simulation after replication ",
    resume$completed_replications, "/", n_replications, "."
  )
} else {
  message("Starting ", n_replications, " LGCM predictor-type replications.")
}
message("Checkpoint file: ", normalizePath(results_file, mustWork = FALSE))

replications_to_run <- if (resume$next_replication <= n_replications) {
  seq.int(resume$next_replication, n_replications)
} else {
  integer()
}
for (replication in replications_to_run) {
  replication_results[[replication]] <- tryCatch(
    run_replication(replication),
    error = function(error) mc_failure(
      replication, "replication", conditionMessage(error)
    )
  )
  progress_counts <- mc_report_progress(
    replication_results[[replication]], replication, n_replications,
    progress_counts
  )
  if (replication %% checkpoint_every == 0L || replication == n_replications) {
    openmx_lgcm_predictor_type_mc <- mc_write_checkpoint(
      replication_results, settings, results_file,
      elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
    )
    message("Checkpoint saved after replication ", replication, ".")
  }
}

openmx_lgcm_predictor_type_mc <- mc_write_checkpoint(
  replication_results, settings, results_file,
  elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
)
message("Simulation complete in ", format(openmx_lgcm_predictor_type_mc$elapsed), ".")
print(openmx_lgcm_predictor_type_mc$summary, row.names = FALSE)
print(openmx_lgcm_predictor_type_mc$test_summary, row.names = FALSE)
print(openmx_lgcm_predictor_type_mc$diagnostic_summary, row.names = FALSE)
