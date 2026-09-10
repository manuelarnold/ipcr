# MGSEM comparison for standard and iterated IPCR in the discrete-time
# autoregressive panel model from Arnold et al. (2020), Simulation I.
#
# The only IPC predictor is a balanced 0/1 group indicator. All ten model
# parameters differ between the two populations in Table 1. The oracle MGSEM
# therefore estimates those ten parameters group-wise while retaining their
# equality constraints across waves. Its joint 10-df likelihood-ratio test is
# compared with the order-invariant Pillai MANCOVA test of the group predictor
# from standard and iterated IPCR.
#
# Run from the package root. For a short smoke run:
# IPCR_MC_REPLICATIONS=2 Rscript tests/manual/Monte-Carlo-Simulation/Simulation-1/lavaan-ar-panel-mgsem-monte-carlo.R


# Settings and package loading --------------------------------------------

group_size <- 125L
sample_size <- 2L * group_size
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
    "lavaan-ar-panel-mgsem-n250-results.rds"
  )
)

if (is.na(n_replications) || n_replications < 1L ||
    is.na(checkpoint_every) || checkpoint_every < 1L) {
  stop("Replication and checkpoint counts must be positive integers.")
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
  group <- sample(rep(0:1, each = group_size))
  model_data <- mc_ar_simulate_data(group, n_waves)
  predictors <- data.frame(group = group)
  multigroup_data <- model_data
  multigroup_data$group <- group

  pooled_attempt <- mc_capture({
    fit <- mc_ar_fit_pooled(model_data, n_waves)
    if (!mc_ar_lavaan_converged(fit)) {
      stop("The pooled lavaan model did not converge.")
    }
    mc_ar_validate_pooled_parameters(fit)
    fit
  })
  if (inherits(pooled_attempt$value, "error")) {
    return(mc_failure(
      replication, "pooled_model", conditionMessage(pooled_attempt$value)
    ))
  }
  pooled_fit <- pooled_attempt$value
  estimates <- diagnostics <- tests <- failures <- list()

  for (method in c("standard", "iterated")) {
    analysis_attempt <- if (method == "standard") {
      mc_capture(suppressMessages(ipcr::ipcr(
        pooled_fit, predictors = predictors, analytic = TRUE
      )))
    } else {
      mc_capture(suppressMessages(ipcr::ipcr_it(
        pooled_fit,
        predictors = predictors,
        analytic = TRUE,
        iteration_info = FALSE
      )))
    }
    if (inherits(analysis_attempt$value, "error")) {
      failures[[length(failures) + 1L]] <- mc_failure(
        replication, "ipcr_analysis", conditionMessage(analysis_attempt$value),
        condition = "dummy_only", method = method,
        error = analysis_attempt$value
      )$failures
      next
    }

    result <- analysis_attempt$value
    estimates[[length(estimates) + 1L]] <- mc_extract_ipcr_estimates(
      result, replication, "dummy_only", method,
      mc_ar_target_intercepts, mc_ar_target_effects,
      informative_name = "group"
    )
    pillai <- mc_pillai_test(result, "group")
    tests[[length(tests) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = method,
      test = "global_heterogeneity",
      statistic = pillai[["statistic"]],
      df1 = pillai[["df1"]],
      df2 = pillai[["df2"]],
      p_value = pillai[["p_value"]],
      row.names = NULL
    )
    is_iterated <- method == "iterated"
    iteration_diagnostics <- mc_iteration_diagnostics(result, is_iterated)
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = method,
      pooled_fit_seconds = pooled_attempt$elapsed_seconds,
      elapsed_seconds = analysis_attempt$elapsed_seconds,
      converged = iteration_diagnostics$converged,
      learning_rate = iteration_diagnostics$learning_rate,
      iterations = iteration_diagnostics$iterations,
      total_iterations = iteration_diagnostics$total_iterations,
      learning_rate_attempts = iteration_diagnostics$learning_rate_attempts,
      C_SE = iteration_diagnostics$C_SE,
      C_l = iteration_diagnostics$C_l,
      attempts = I(list(iteration_diagnostics$attempts)),
      pooled_fit_warnings = paste(pooled_attempt$warnings, collapse = " | "),
      warnings = paste(analysis_attempt$warnings, collapse = " | "),
      row.names = NULL
    )
  }

  full_attempt <- mc_capture({
    fit <- mc_ar_fit_multigroup(
      multigroup_data, constrained = FALSE, n_waves = n_waves
    )
    if (!mc_ar_lavaan_converged(fit)) {
      stop("The unrestricted multigroup lavaan model did not converge.")
    }
    fit
  })
  null_attempt <- mc_capture({
    fit <- mc_ar_fit_multigroup(
      multigroup_data, constrained = TRUE, n_waves = n_waves
    )
    if (!mc_ar_lavaan_converged(fit)) {
      stop("The equality-constrained multigroup lavaan model did not converge.")
    }
    fit
  })
  if (inherits(full_attempt$value, "error") ||
      inherits(null_attempt$value, "error")) {
    failed_attempt <- if (inherits(full_attempt$value, "error")) {
      full_attempt$value
    } else null_attempt$value
    failures[[length(failures) + 1L]] <- mc_failure(
      replication, "mgsem", conditionMessage(failed_attempt),
      condition = "dummy_only", method = "MGSEM"
    )$failures
  } else {
    full_fit <- full_attempt$value
    null_fit <- null_attempt$value
    estimates[[length(estimates) + 1L]] <- mc_ar_extract_mgsem_estimates(
      full_fit, replication
    )
    tests[[length(tests) + 1L]] <- mc_ar_extract_lrt(
      null_fit, full_fit, replication
    )
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = "MGSEM",
      pooled_fit_seconds = pooled_attempt$elapsed_seconds,
      elapsed_seconds = full_attempt$elapsed_seconds +
        null_attempt$elapsed_seconds,
      converged = NA,
      learning_rate = NA_real_,
      iterations = NA_integer_,
      total_iterations = NA_integer_,
      learning_rate_attempts = NA_integer_,
      C_SE = NA_real_,
      C_l = NA_real_,
      attempts = I(list(data.frame())),
      pooled_fit_warnings = paste(pooled_attempt$warnings, collapse = " | "),
      warnings = paste(
        c(full_attempt$warnings, null_attempt$warnings), collapse = " | "
      ),
      row.names = NULL
    )
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
    "Arnold et al. (2020) discrete-time autoregressive panel:",
    "IPCR versus oracle MGSEM"
  ),
  sample_size = sample_size,
  group_size = group_size,
  n_waves = n_waves,
  n_replications = n_replications,
  checkpoint_every = checkpoint_every,
  seed = simulation_seed,
  conditions = "dummy_only",
  methods = c("standard", "iterated", "MGSEM"),
  ipc_predictor = "balanced group dummy only",
  mgsem_group_order = c("0", "1"),
  mgsem_group_specific_parameters = names(mc_ar_target_effects),
  mgsem_equal_across_waves = TRUE,
  joint_test = paste(
    "IPCR Pillai MANCOVA versus 10-parameter MGSEM likelihood-ratio test"
  ),
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
    "Resuming the AR-panel MGSEM-comparison simulation after replication ",
    resume$completed_replications, "/", n_replications, "."
  )
} else {
  message("Starting ", n_replications, " AR-panel MGSEM-comparison replications.")
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
    lavaan_ar_panel_mgsem_mc <- mc_write_checkpoint(
      replication_results, settings, results_file,
      elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
    )
    message("Checkpoint saved after replication ", replication, ".")
  }
}

lavaan_ar_panel_mgsem_mc <- mc_write_checkpoint(
  replication_results, settings, results_file,
  elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
)
message("Simulation complete in ", format(lavaan_ar_panel_mgsem_mc$elapsed), ".")
print(lavaan_ar_panel_mgsem_mc$summary, row.names = FALSE)
print(lavaan_ar_panel_mgsem_mc$test_summary, row.names = FALSE)
print(lavaan_ar_panel_mgsem_mc$diagnostic_summary, row.names = FALSE)
