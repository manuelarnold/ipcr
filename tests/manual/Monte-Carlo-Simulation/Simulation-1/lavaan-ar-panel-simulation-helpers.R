# Shared population model, data generator, and lavaan model builders for the
# two autoregressive-panel Monte Carlo scripts. This file is sourced by the
# runnable simulations; it does not run a simulation itself.

mc_ar_group_0_parameters <- list(
  transition = matrix(
    c(0.700, 0.000,
      0.000, 0.700),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  ),
  initial_covariance = matrix(
    c(1.000, 0.300,
      0.300, 1.000),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  ),
  innovation_covariance = matrix(
    c(0.510, 0.153,
      0.153, 0.510),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  )
)

mc_ar_group_1_parameters <- list(
  transition = matrix(
    c(0.450, 0.300,
      0.300, 0.450),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  ),
  initial_covariance = matrix(
    c(2.000, 1.000,
      1.000, 2.000),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  ),
  innovation_covariance = matrix(
    c(1.145, 0.168,
      0.168, 1.145),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("x", "y"))
  )
)

mc_ar_target_intercepts <- c(
  beta_xx = 0.700,
  beta_xy = 0.000,
  beta_yx = 0.000,
  beta_yy = 0.700,
  phi_xx = 1.000,
  phi_yx = 0.300,
  phi_yy = 1.000,
  psi_xx = 0.510,
  psi_yx = 0.153,
  psi_yy = 0.510
)

mc_ar_target_effects <- c(
  beta_xx = -0.250,
  beta_xy = 0.300,
  beta_yx = 0.300,
  beta_yy = -0.250,
  phi_xx = 1.000,
  phi_yx = 0.700,
  phi_yy = 1.000,
  psi_xx = 0.635,
  psi_yx = 0.015,
  psi_yy = 0.635
)

mc_ar_spectral_radius <- function(transition) {
  max(Mod(eigen(transition, only.values = TRUE)$values))
}

mc_ar_validate_population <- function() {
  interpolation_grid <- seq(0, 1, length.out = 101L)
  radii <- vapply(interpolation_grid, function(informative) {
    transition <- mc_ar_group_0_parameters$transition + informative *
      (mc_ar_group_1_parameters$transition -
         mc_ar_group_0_parameters$transition)
    mc_ar_spectral_radius(transition)
  }, numeric(1))
  if (any(radii >= 1)) {
    stop("The interpolated transition model is not dynamically stable.")
  }
  invisible(TRUE)
}

# Draw zero-mean bivariate normal observations with a covariance matrix that
# may differ by row. This avoids repeatedly decomposing 2-by-2 matrices.
mc_ar_rmvnorm_varying <- function(variance_x, covariance_xy, variance_y) {
  n <- length(variance_x)
  standard_x <- stats::rnorm(n)
  standard_y <- stats::rnorm(n)
  sd_x <- sqrt(variance_x)
  conditional_variance <- variance_y - covariance_xy^2 / variance_x
  if (any(!is.finite(conditional_variance)) ||
      any(conditional_variance <= 0)) {
    stop("An interpolated covariance matrix is not positive definite.")
  }
  cbind(
    x = sd_x * standard_x,
    y = covariance_xy / sd_x * standard_x +
      sqrt(conditional_variance) * standard_y
  )
}

mc_ar_interpolate <- function(informative, group_0, group_1) {
  group_0 + informative * (group_1 - group_0)
}

# The informative predictor lies in [0, 1]. At the endpoints this reproduces
# the two groups in Table 1; intermediate values linearly interpolate every
# heterogeneous model parameter.
mc_ar_simulate_data <- function(informative, n_waves = 5L) {
  if (any(!is.finite(informative)) ||
      any(informative < 0 | informative > 1)) {
    stop("The informative predictor must contain values in [0, 1].")
  }
  n <- length(informative)
  interpolate_entry <- function(component, row, column) {
    mc_ar_interpolate(
      informative,
      mc_ar_group_0_parameters[[component]][row, column],
      mc_ar_group_1_parameters[[component]][row, column]
    )
  }

  state <- mc_ar_rmvnorm_varying(
    interpolate_entry("initial_covariance", 1L, 1L),
    interpolate_entry("initial_covariance", 1L, 2L),
    interpolate_entry("initial_covariance", 2L, 2L)
  )
  observations <- matrix(
    NA_real_,
    nrow = n,
    ncol = 2L * n_waves,
    dimnames = list(
      NULL,
      as.vector(rbind(
        paste0("x", seq_len(n_waves)),
        paste0("y", seq_len(n_waves))
      ))
    )
  )
  observations[, c("x1", "y1")] <- state

  beta_xx <- interpolate_entry("transition", 1L, 1L)
  beta_xy <- interpolate_entry("transition", 1L, 2L)
  beta_yx <- interpolate_entry("transition", 2L, 1L)
  beta_yy <- interpolate_entry("transition", 2L, 2L)
  psi_xx <- interpolate_entry("innovation_covariance", 1L, 1L)
  psi_yx <- interpolate_entry("innovation_covariance", 1L, 2L)
  psi_yy <- interpolate_entry("innovation_covariance", 2L, 2L)

  for (wave in 2:n_waves) {
    innovations <- mc_ar_rmvnorm_varying(psi_xx, psi_yx, psi_yy)
    previous_x <- state[, "x"]
    previous_y <- state[, "y"]
    state <- cbind(
      x = beta_xx * previous_x + beta_xy * previous_y + innovations[, "x"],
      y = beta_yx * previous_x + beta_yy * previous_y + innovations[, "y"]
    )
    observations[, c(paste0("x", wave), paste0("y", wave))] <- state
  }
  as.data.frame(observations)
}

mc_ar_lavaan_label <- function(parameter, multigroup, constrained) {
  if (!multigroup) return(parameter)
  if (constrained) {
    paste0("c(", parameter, ",", parameter, ")")
  } else {
    paste0("c(", parameter, "_g0,", parameter, "_g1)")
  }
}

mc_ar_make_lavaan_syntax <- function(
    n_waves = 5L, multigroup = FALSE, constrained = FALSE) {
  label <- function(parameter) {
    mc_ar_lavaan_label(parameter, multigroup, constrained)
  }
  transitions <- unlist(lapply(2:n_waves, function(wave) {
    previous <- wave - 1L
    c(
      paste0(
        "x", wave, " ~ ", label("beta_xx"), "*x", previous,
        " + ", label("beta_xy"), "*y", previous
      ),
      paste0(
        "y", wave, " ~ ", label("beta_yx"), "*x", previous,
        " + ", label("beta_yy"), "*y", previous
      )
    )
  }))
  innovation_covariances <- unlist(lapply(2:n_waves, function(wave) {
    c(
      paste0("x", wave, " ~~ ", label("psi_xx"), "*x", wave),
      paste0("x", wave, " ~~ ", label("psi_yx"), "*y", wave),
      paste0("y", wave, " ~~ ", label("psi_yy"), "*y", wave)
    )
  }))
  paste(
    c(
      transitions,
      paste0("x1 ~~ ", label("phi_xx"), "*x1"),
      paste0("x1 ~~ ", label("phi_yx"), "*y1"),
      paste0("y1 ~~ ", label("phi_yy"), "*y1"),
      innovation_covariances
    ),
    collapse = "\n"
  )
}

mc_ar_fit_pooled <- function(data, n_waves = 5L) {
  lavaan::sem(
    model = mc_ar_make_lavaan_syntax(n_waves),
    data = data,
    estimator = "ML",
    meanstructure = FALSE,
    fixed.x = FALSE,
    warn = TRUE
  )
}

mc_ar_fit_multigroup <- function(data, constrained, n_waves = 5L) {
  lavaan::sem(
    model = mc_ar_make_lavaan_syntax(
      n_waves,
      multigroup = TRUE,
      constrained = constrained
    ),
    data = data,
    group = "group",
    # Do not let the first observed row determine which population receives
    # the _g0 and _g1 parameter labels.
    group.label = c("0", "1"),
    estimator = "ML",
    meanstructure = FALSE,
    fixed.x = FALSE,
    warn = TRUE
  )
}

mc_ar_lavaan_converged <- function(fit) {
  isTRUE(lavaan::lavInspect(fit, "converged"))
}

mc_ar_validate_pooled_parameters <- function(fit) {
  fitted_names <- names(ipcr:::coef_ipcr(fit))
  expected_names <- names(mc_ar_target_effects)
  if (!identical(fitted_names, expected_names)) {
    stop(
      "Unexpected reduced lavaan parameter names or order: ",
      paste(fitted_names, collapse = ", ")
    )
  }
  invisible(TRUE)
}

# lavaan retains duplicated free coordinates when equality is imposed through
# labels. Reduce estimates and their covariance matrix to one coordinate for
# each unique label, just as ipcr does for its score and bread calculations.
mc_ar_reduced_lavaan_estimates <- function(fit) {
  estimates <- lavaan::coef(fit)
  covariance <- lavaan::vcov(fit)
  reduction <- ipcr:::lavaan_constraint_reduction_ipcr(fit@Model)
  if (!is.null(reduction)) {
    estimates <- estimates[reduction$keep]
    covariance <- covariance[
      reduction$keep, reduction$keep, drop = FALSE
    ]
  }
  if (anyDuplicated(names(estimates))) {
    stop("The reduced lavaan coefficient names are not unique.")
  }
  dimnames(covariance) <- list(names(estimates), names(estimates))
  list(estimates = estimates, covariance = covariance)
}

mc_ar_extract_mgsem_estimates <- function(fit, replication) {
  group_labels <- as.character(lavaan::lavInspect(fit, "group.label"))
  if (!identical(group_labels, c("0", "1"))) {
    stop("The lavaan group order is not the expected 0, 1 order.")
  }
  reduced <- mc_ar_reduced_lavaan_estimates(fit)
  estimates <- reduced$estimates
  covariance <- reduced$covariance
  expected_labels <- unlist(lapply(
    names(mc_ar_target_effects),
    function(parameter) paste0(parameter, c("_g0", "_g1"))
  ))
  if (!setequal(names(estimates), expected_labels)) {
    stop("The unrestricted MGSEM does not contain the expected 20 parameters.")
  }
  critical_value <- stats::qnorm(0.975)

  rows <- lapply(names(mc_ar_target_effects), function(parameter) {
    labels <- paste0(parameter, c("_g0", "_g1"))
    if (!all(labels %in% names(estimates))) {
      stop("Unable to find MGSEM parameters for ", parameter, ".")
    }
    contrast <- c(-1, 1)
    contrast_covariance <- covariance[labels, labels, drop = FALSE]
    group_0_estimate <- unname(estimates[labels[1L]])
    group_0_se <- sqrt(covariance[labels[1L], labels[1L]])
    group_effect <- sum(contrast * estimates[labels])
    group_effect_se <- sqrt(as.numeric(
      t(contrast) %*% contrast_covariance %*% contrast
    ))
    row_estimates <- c(group_0_estimate, group_effect)
    row_standard_errors <- c(group_0_se, group_effect_se)
    data.frame(
      replication = replication,
      condition = "dummy_only",
      method = "MGSEM",
      parameter = parameter,
      predictor = c("(Intercept)", "group"),
      target = c(
        unname(mc_ar_target_intercepts[[parameter]]),
        unname(mc_ar_target_effects[[parameter]])
      ),
      estimate = row_estimates,
      standard_error = row_standard_errors,
      residual_df = NA_real_,
      p_value = 2 * stats::pnorm(
        -abs(row_estimates / row_standard_errors)
      ),
      ci_lower = row_estimates - critical_value * row_standard_errors,
      ci_upper = row_estimates + critical_value * row_standard_errors,
      row.names = NULL
    )
  })
  mc_combine_rows(rows)
}

mc_ar_extract_lrt <- function(null_fit, full_fit, replication) {
  comparison <- suppressMessages(lavaan::lavTestLRT(null_fit, full_fit))
  comparison <- as.data.frame(comparison)
  statistic_column <- intersect(c("Chisq diff", "Chisq.diff"), names(comparison))
  df_column <- intersect(c("Df diff", "Df.diff"), names(comparison))
  p_column <- grep("^Pr\\(", names(comparison), value = TRUE)
  if (length(statistic_column) != 1L || length(df_column) != 1L ||
      length(p_column) != 1L) {
    stop("Unable to identify the lavaan likelihood-ratio columns.")
  }
  comparison_row <- which(is.finite(comparison[[statistic_column]]))
  if (length(comparison_row) != 1L) {
    stop("Unable to identify the lavaan likelihood-ratio comparison row.")
  }
  comparison_df <- comparison[[df_column]][comparison_row]
  if (!identical(as.integer(comparison_df), length(mc_ar_target_effects))) {
    stop("The joint MGSEM likelihood-ratio test does not have 10 df.")
  }
  data.frame(
    replication = replication,
    condition = "dummy_only",
    method = "MGSEM",
    test = "global_heterogeneity",
    statistic = comparison[[statistic_column]][comparison_row],
    df1 = comparison_df,
    df2 = NA_real_,
    p_value = comparison[[p_column]][comparison_row],
    row.names = NULL
  )
}
