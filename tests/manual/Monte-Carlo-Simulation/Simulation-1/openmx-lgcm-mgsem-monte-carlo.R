# MGSEM comparison for standard and iterated IPCR in an OpenMx latent growth
# curve model (LGCM).
#
# The only IPC predictor is a balanced 0/1 group indicator. The groups differ
# only in the manifest residual variance (0.3 versus 0.5). The oracle MGSEM
# estimates that residual variance separately while constraining all other
# parameters equal. Its one-df likelihood-ratio test is compared with both the
# IPCR t-test of the residual-variance slope and the global Pillai MANCOVA.
#
# Run from the package root. For a short smoke run:
# IPCR_MC_REPLICATIONS=2 Rscript tests/manual/Monte-Carlo-Simulation/Simulation-1/openmx-lgcm-mgsem-monte-carlo.R


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
  file.path(simulation_directory, "Results", "openmx-lgcm-mgsem-n250-results.rds")
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


# LGCM data generation and model construction -----------------------------

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

simulate_lgcm_data <- function(group) {
  data <- matrix(
    NA_real_, nrow = length(group), ncol = length(time_scores),
    dimnames = list(NULL, manifest_variables)
  )
  residual_variance <- 0.3 + 0.2 * group
  for (individual in seq_along(group)) {
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

make_lgcm_model <- function(data, model_name, residual_label) {
  OpenMx::mxModel(
    model_name,
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
      values = 0.4, labels = residual_label, lbound = 0.0001
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

make_multigroup_model <- function(data_group_0, data_group_1, constrained) {
  residual_labels <- if (constrained) {
    c("residual", "residual")
  } else {
    c("residual_g0", "residual_g1")
  }
  group_0_model <- make_lgcm_model(
    data_group_0, "group0", residual_labels[1L]
  )
  group_1_model <- make_lgcm_model(
    data_group_1, "group1", residual_labels[2L]
  )
  OpenMx::mxModel(
    if (constrained) "MGSEMNull" else "MGSEMFull",
    group_0_model,
    group_1_model,
    OpenMx::mxFitFunctionMultigroup(
      c("group0.fitfunction", "group1.fitfunction")
    )
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

extract_openmx_mgsem <- function(fit, replication) {
  estimates <- OpenMx::omxGetParameters(fit)
  covariance <- fit$output$vcov
  required <- c("residual_g0", "residual_g1")
  if (!all(required %in% names(estimates)) ||
      !all(required %in% rownames(covariance))) {
    stop("Unable to extract the OpenMx MGSEM residual contrast.")
  }
  contrast <- c(-1, 1)
  contrast_covariance <- covariance[required, required, drop = FALSE]
  group_effect <- sum(contrast * estimates[required])
  group_effect_se <- sqrt(as.numeric(
    t(contrast) %*% contrast_covariance %*% contrast
  ))
  group_0_se <- sqrt(covariance["residual_g0", "residual_g0"])
  critical_value <- stats::qnorm(0.975)

  data.frame(
    replication = replication,
    condition = "dummy_only",
    method = "MGSEM",
    parameter = "residual",
    predictor = c("(Intercept)", "group"),
    target = c(0.3, 0.2),
    estimate = c(estimates[["residual_g0"]], group_effect),
    standard_error = c(group_0_se, group_effect_se),
    residual_df = NA_real_,
    p_value = 2 * stats::pnorm(
      -abs(c(estimates[["residual_g0"]] / group_0_se,
             group_effect / group_effect_se))
    ),
    ci_lower = c(estimates[["residual_g0"]], group_effect) -
      critical_value * c(group_0_se, group_effect_se),
    ci_upper = c(estimates[["residual_g0"]], group_effect) +
      critical_value * c(group_0_se, group_effect_se),
    row.names = NULL
  )
}


# One replication ---------------------------------------------------------

run_replication <- function(replication) {
  set.seed(simulation_seed + replication - 1L)
  group <- rep(0:1, each = group_size)
  model_data <- simulate_lgcm_data(group)
  predictors <- data.frame(group = group)
  data_group_0 <- model_data[group == 0, , drop = FALSE]
  data_group_1 <- model_data[group == 1, , drop = FALSE]

  pooled_attempt <- mc_capture(fit_openmx(
    make_lgcm_model(model_data, "PooledLGCM", "residual")
  ))
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
      target_intercepts, target_effects, informative_name = "group"
    )
    residual_test <- estimates[[length(estimates)]]
    residual_test <- residual_test[
      residual_test$parameter == "residual" &
        residual_test$predictor == "group",
      ,
      drop = FALSE
    ]
    tests[[length(tests) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = method,
      test = "residual_group_effect",
      statistic = residual_test$estimate / residual_test$standard_error,
      df1 = 1,
      df2 = residual_test$residual_df,
      p_value = residual_test$p_value,
      row.names = NULL
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

  full_attempt <- mc_capture(fit_openmx(make_multigroup_model(
    data_group_0, data_group_1, constrained = FALSE
  )))
  null_attempt <- mc_capture(fit_openmx(make_multigroup_model(
    data_group_0, data_group_1, constrained = TRUE
  )))
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
    estimates[[length(estimates) + 1L]] <- extract_openmx_mgsem(
      full_fit, replication
    )
    comparison <- suppressMessages(OpenMx::mxCompare(
      base = full_fit,
      comparison = null_fit
    ))
    comparison_row <- comparison[2L, , drop = FALSE]
    tests[[length(tests) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = "MGSEM",
      test = "global_heterogeneity",
      statistic = comparison_row$diffLL,
      df1 = comparison_row$diffdf,
      df2 = NA_real_,
      p_value = comparison_row$p,
      row.names = NULL
    )
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      replication = replication,
      condition = "dummy_only",
      method = "MGSEM",
      pooled_fit_seconds = pooled_attempt$elapsed_seconds,
      elapsed_seconds = full_attempt$elapsed_seconds + null_attempt$elapsed_seconds,
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
        c(full_attempt$warnings, null_attempt$warnings),
        collapse = " | "
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
  design = "OpenMx LGCM IPCR versus oracle MGSEM",
  sample_size = sample_size,
  group_size = group_size,
  n_replications = n_replications,
  checkpoint_every = checkpoint_every,
  seed = simulation_seed,
  conditions = "dummy_only",
  methods = c("standard", "iterated", "MGSEM"),
  ipc_predictor = "balanced group dummy only",
  mgsem_group_specific_parameters = "residual",
  mgsem_equal_parameters = setdiff(names(target_effects), "residual"),
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
    "Resuming the LGCM MGSEM-comparison simulation after replication ",
    resume$completed_replications, "/", n_replications, "."
  )
} else {
  message("Starting ", n_replications, " LGCM MGSEM-comparison replications.")
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
    openmx_lgcm_mgsem_mc <- mc_write_checkpoint(
      replication_results, settings, results_file,
      elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
    )
    message("Checkpoint saved after replication ", replication, ".")
  }
}

openmx_lgcm_mgsem_mc <- mc_write_checkpoint(
  replication_results, settings, results_file,
  elapsed = mc_elapsed_since(started_at, resume$elapsed_seconds)
)
message("Simulation complete in ", format(openmx_lgcm_mgsem_mc$elapsed), ".")
print(openmx_lgcm_mgsem_mc$summary, row.names = FALSE)
print(openmx_lgcm_mgsem_mc$test_summary, row.names = FALSE)
print(openmx_lgcm_mgsem_mc$diagnostic_summary, row.names = FALSE)
