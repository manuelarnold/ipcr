# Predictor-type Monte Carlo simulation for the discrete-time autoregressive
# panel model from Arnold et al. (2020), Simulation I.
#
# The informative predictor is either a balanced 0/1 dummy or U(0, 1). Every
# model also includes two independent noise predictors: one balanced dummy and
# one U(0, 1) variable. Values between zero and one linearly interpolate all
# ten heterogeneous parameters between the two populations in Table 1. The
# transition matrix remains stationary throughout this interval.
#
# Standard and iterated IPCR use analytical lavaan scores. Run from the package
# root. For a short smoke run:
# IPCR_MC_REPLICATIONS=2 Rscript tests/manual/Monte-Carlo-Simulation/Simulation-1/lavaan-ar-panel-predictor-type-monte-carlo.R


# Settings and package loading --------------------------------------------

sample_size <- 250L
group_size <- sample_size %/% 2L
n_waves <- 5L
n_replications <- as.integer(Sys.getenv("IPCR_MC_REPLICATIONS", "1000"))
checkpoint_every <- as.integer(Sys.getenv("IPCR_MC_CHECKPOINT_EVERY", "10"))
simulation_seed <- 6202020L
simulation_directory <- file.path(
  "tests", "manual", "Monte-Carlo-Simulation", "Simulation-1"
)
results_file <- Sys.getenv(
  "IPCR_MC_RESULTS_FILE",
  file.path(
    simulation_directory, "Results",
    "lavaan-ar-panel-predictor-type-n250-results.rds"
  )
)

if (sample_size %% 2L != 0L || is.na(n_replications) ||
    n_replications < 1L || is.na(checkpoint_every) ||
    checkpoint_every < 1L) {
  stop("The sample size must be even and counts must be positive integers.")
}
required_packages <- c("lavaan", "pkgload")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
)]
if (length(missing_packages) > 0L) {
  stop("This simulation requires: ", paste(missing_packages, collapse = ", "))
}
if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] != "ipcr") {
  stop("Run this script from the ipcr package root.")
}

pkgload::load_all(".", quiet = TRUE)
iterated_convergence <- eval(formals(ipcr_it)$conv)
iterated_learning_rate <- eval(formals(ipcr_it)$learning_rate)
iterated_max_iterations <- eval(formals(ipcr_it)$max_it)
source(file.path(simulation_directory, "ipcr-monte-carlo-helpers.R"))
source(file.path(simulation_directory, "lavaan-ar-panel-simulation-helpers.R"))
mc_ar_validate_population()


# One replication ---------------------------------------------------------

run_replication <- function(replication) {
  set.seed(simulation_seed + replication - 1L)
  conditions <- list(
    informative_dummy = sample(rep(0:1, each = group_size)),
    informative_continuous = stats::runif(sample_size, 0, 1)
  )
  estimates <- diagnostics <- tests <- failures <- list()

  for (condition in names(conditions)) {
    informative <- conditions[[condition]]
    predictors <- data.frame(
      informative = informative,
      noise_dummy = sample(rep(0:1, each = group_size)),
      noise_uniform = stats::runif(sample_size, 0, 1)
    )
    model_data <- mc_ar_simulate_data(informative, n_waves)
    fit_attempt <- mc_capture({
      fit <- mc_ar_fit_pooled(model_data, n_waves)
      if (!mc_ar_lavaan_converged(fit)) {
        stop("The pooled lavaan model did not converge.")
      }
      mc_ar_validate_pooled_parameters(fit)
      fit
    })
    if (inherits(fit_attempt$value, "error")) {
      failures[[length(failures) + 1L]] <- mc_failure(
        replication, "pooled_model", conditionMessage(fit_attempt$value),
        condition = condition
      )$failures
      next
    }
    fit <- fit_attempt$value

    for (method in c("standard", "iterated")) {
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
      estimates[[length(estimates) + 1L]] <- mc_extract_ipcr_estimates(
        result, replication, condition, method,
        mc_ar_target_intercepts, mc_ar_target_effects,
        informative_name = "informative"
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
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
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
  design = paste(
    "Arnold et al. (2020) discrete-time autoregressive panel",
    "predictor-type simulation"
  ),
  sample_size = sample_size,
  n_waves = n_waves,
  n_replications = n_replications,
  checkpoint_every = checkpoint_every,
  seed = simulation_seed,
  conditions = c("informative_dummy", "informative_continuous"),
  methods = c("standard", "iterated"),
  informative_continuous_distribution = "Uniform(0, 1)",
  noise_predictors = c("balanced dummy", "Uniform(0, 1)"),
  parameter_interpolation = "group 0 + informative * (group 1 - group 0)",
  group_0_parameters = mc_ar_group_0_parameters,
  group_1_parameters = mc_ar_group_1_parameters,
  target_intercepts = mc_ar_target_intercepts,
  target_effects = mc_ar_target_effects,
  analytic = TRUE,
  iterated_convergence = iterated_convergence,
  iterated_learning_rate = iterated_learning_rate,
  iterated_max_iterations = iterated_max_iterations,
  ipcr_version = as.character(utils::packageVersion("ipcr")),
  lavaan_version = as.character(utils::packageVersion("lavaan"))
)

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
resume <- mc_resume_checkpoint(results_file, settings, n_replications)
replication_results <- resume$results
progress_counts <- resume$counts
started_at <- Sys.time()
if (resume$completed_replications > 0L) {
  message(
    "Resuming the AR-panel predictor-type simulation after replication ",
    resume$completed_replications, "/", n_replications, "."
  )
} else {
  message("Starting ", n_replications, " AR-panel predictor-type replications.")
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
    lavaan_ar_panel_predictor_type_mc <- mc_write_checkpoint(
      replication_results, settings, results_file,
      elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
    )
    message("Checkpoint saved after replication ", replication, ".")
  }
}

lavaan_ar_panel_predictor_type_mc <- mc_write_checkpoint(
  replication_results, settings, results_file,
  elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
)
message(
  "Simulation complete in ",
  format(lavaan_ar_panel_predictor_type_mc$elapsed), "."
)
print(lavaan_ar_panel_predictor_type_mc$summary, row.names = FALSE)
print(lavaan_ar_panel_predictor_type_mc$test_summary, row.names = FALSE)
print(lavaan_ar_panel_predictor_type_mc$diagnostic_summary, row.names = FALSE)
