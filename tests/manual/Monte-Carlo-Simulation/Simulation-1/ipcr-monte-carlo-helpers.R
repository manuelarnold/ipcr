# Shared infrastructure for the manual IPCR Monte Carlo scripts.
#
# This file defines extraction and checkpoint helpers only. It does not run a
# simulation by itself. Simulation scripts source it from the package root.

mc_empty_data_frame <- function() data.frame()

mc_combine_rows <- function(x) {
  x <- Filter(function(element) !is.null(element) && NROW(element) > 0L, x)
  if (length(x) == 0L) return(mc_empty_data_frame())

  # Checkpoints written by earlier script versions can lack newer diagnostic
  # columns. Fill those columns before row-binding so a compatible partial run
  # can still be resumed.
  column_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  list_columns <- setNames(
    vapply(column_names, function(column_name) {
      any(vapply(x, function(element) {
        column_name %in% names(element) && is.list(element[[column_name]])
      }, logical(1)))
    }, logical(1)),
    column_names
  )
  x <- lapply(x, function(element) {
    missing_columns <- setdiff(column_names, names(element))
    for (column_name in missing_columns) {
      element[[column_name]] <- if (list_columns[[column_name]]) {
        I(rep(list(NULL), NROW(element)))
      } else {
        rep(NA, NROW(element))
      }
    }
    element[column_names]
  })
  do.call(rbind, x)
}

mc_capture <- function(expression) {
  captured_warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(warning) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(warning)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  list(
    value = value,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
    warnings = unique(captured_warnings)
  )
}

mc_iteration_diagnostics <- function(result, is_iterated) {
  empty_attempts <- data.frame(
    learning_rate = numeric(),
    iterations = integer(),
    C_SE = numeric(),
    C_l = numeric(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(is_iterated)) {
    return(list(
      converged = NA,
      learning_rate = NA_real_,
      iterations = NA_integer_,
      total_iterations = NA_integer_,
      learning_rate_attempts = NA_integer_,
      C_SE = NA_real_,
      C_l = NA_real_,
      attempts = empty_attempts
    ))
  }
  attempts <- result$convergence$attempts
  if (is.null(attempts)) attempts <- empty_attempts
  scalar_or_na <- function(value, mode = c("numeric", "integer")) {
    mode <- match.arg(mode)
    if (length(value) != 1L || is.na(value)) {
      return(if (mode == "integer") NA_integer_ else NA_real_)
    }
    if (mode == "integer") as.integer(value) else as.numeric(value)
  }
  list(
    converged = isTRUE(result$convergence$converged),
    learning_rate = scalar_or_na(result$convergence$learning_rate),
    iterations = scalar_or_na(result$info$iterations, "integer"),
    total_iterations = if (NROW(attempts) == 0L ||
        !"iterations" %in% names(attempts)) {
      scalar_or_na(result$info$iterations, "integer")
    } else {
      as.integer(sum(attempts$iterations, na.rm = TRUE))
    },
    learning_rate_attempts = NROW(attempts),
    C_SE = scalar_or_na(result$convergence$C_SE),
    C_l = scalar_or_na(result$convergence$C_l),
    attempts = attempts
  )
}

mc_replication_counts <- function(result) {
  nonconvergence_errors <- sum(
    mc_ipcr_it_nonconvergence_flags(result$failures)
  )
  c(
    failures = NROW(result$failures),
    nonconvergence_errors = nonconvergence_errors
  )
}

mc_report_progress <- function(
    result, replication, n_replications, cumulative_counts) {
  cumulative_counts <- cumulative_counts + mc_replication_counts(result)
  message(
    "Completed ", replication, "/", n_replications,
    " (", cumulative_counts[["failures"]], " failures, including ",
    cumulative_counts[["nonconvergence_errors"]],
    " iterated-IPCR nonconvergence errors)."
  )
  cumulative_counts
}

mc_extract_ipcr_estimates <- function(
    result, replication, condition, method, target_intercepts,
    target_effects, informative_name) {
  coefficients <- stats::coef(result)
  parameter_names <- colnames(coefficients)
  coefficient_names <- rownames(coefficients)
  regression_summaries <- summary(result$mlm)
  names(regression_summaries) <- sub(
    "^Response ", "", names(regression_summaries)
  )
  residual_df <- stats::df.residual(result$mlm)
  critical_value <- stats::qt(0.975, df = residual_df)

  rows <- lapply(parameter_names, function(parameter) {
    lapply(coefficient_names, function(coefficient) {
      coefficient_table <- regression_summaries[[parameter]]$coefficients
      estimate <- coefficient_table[coefficient, "Estimate"]
      standard_error <- coefficient_table[coefficient, "Std. Error"]
      p_value <- coefficient_table[coefficient, "Pr(>|t|)"]
      target <- if (identical(coefficient, "(Intercept)")) {
        unname(target_intercepts[[parameter]])
      } else if (identical(coefficient, informative_name)) {
        unname(target_effects[[parameter]])
      } else {
        0
      }
      data.frame(
        replication = replication,
        condition = condition,
        method = method,
        parameter = parameter,
        predictor = coefficient,
        target = target,
        estimate = unname(estimate),
        standard_error = unname(standard_error),
        residual_df = residual_df,
        p_value = unname(p_value),
        ci_lower = unname(estimate - critical_value * standard_error),
        ci_upper = unname(estimate + critical_value * standard_error),
        row.names = NULL
      )
    })
  })
  mc_combine_rows(unlist(rows, recursive = FALSE))
}

mc_summarize_estimates <- function(estimates) {
  if (NROW(estimates) == 0L) return(mc_empty_data_frame())
  groups <- split(
    estimates,
    interaction(
      estimates$condition,
      estimates$method,
      estimates$parameter,
      estimates$predictor,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  mc_combine_rows(lapply(groups, function(group) {
    valid <- is.finite(group$estimate) &
      is.finite(group$standard_error) &
      is.finite(group$p_value)
    group <- group[valid, , drop = FALSE]
    if (NROW(group) == 0L) return(NULL)
    target <- group$target[1L]
    data.frame(
      condition = group$condition[1L],
      method = group$method[1L],
      parameter = group$parameter[1L],
      predictor = group$predictor[1L],
      target = target,
      mean_estimate = mean(group$estimate),
      bias = mean(group$estimate - target),
      empirical_sd = stats::sd(group$estimate),
      mean_standard_error = mean(group$standard_error),
      rmse = sqrt(mean((group$estimate - target)^2)),
      coverage_95 = mean(
        group$ci_lower <= target & group$ci_upper >= target
      ),
      rejection_rate_05 = mean(group$p_value < 0.05),
      successful_replications = NROW(group),
      row.names = NULL
    )
  }))
}

mc_condition_method_key <- function(condition, method) {
  encode <- function(value) {
    value <- as.character(value)
    ifelse(is.na(value), "<NA>", paste0(nchar(value), ":", value))
  }
  paste(encode(condition), encode(method), sep = "\u001f")
}

mc_result_row_key <- function(x) {
  required <- c("replication", "condition", "method")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Monte Carlo rows are missing fields required for matching: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  paste(
    as.integer(x$replication),
    mc_condition_method_key(x$condition, x$method),
    sep = "\u001e"
  )
}

mc_keep_converged_iterated <- function(x, diagnostics) {
  if (NROW(x) == 0L || !"method" %in% names(x)) return(x)
  is_iterated <- !is.na(x$method) & x$method == "iterated"
  if (!any(is_iterated)) return(x)

  required <- c("replication", "condition", "method", "converged")
  missing <- setdiff(required, names(diagnostics))
  if (length(missing) > 0L) {
    stop(
      "Iterated-IPCR results cannot be matched to convergence diagnostics; ",
      "the diagnostics are missing: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  converged_diagnostics <- diagnostics[
    !is.na(diagnostics$method) & diagnostics$method == "iterated" &
      !is.na(diagnostics$converged) & diagnostics$converged,
    ,
    drop = FALSE
  ]
  converged_keys <- if (NROW(converged_diagnostics) == 0L) {
    character()
  } else {
    unique(mc_result_row_key(converged_diagnostics))
  }
  keep <- !is_iterated
  keep[is_iterated] <- mc_result_row_key(x[is_iterated, , drop = FALSE]) %in%
    converged_keys
  x[keep, , drop = FALSE]
}

mc_collect_attempt_diagnostics <- function(diagnostics) {
  if (NROW(diagnostics) == 0L ||
      !all(c("replication", "condition", "method", "attempts") %in%
        names(diagnostics))) {
    return(mc_empty_data_frame())
  }
  rows <- lapply(seq_len(NROW(diagnostics)), function(index) {
    if (is.na(diagnostics$method[index]) ||
        diagnostics$method[index] != "iterated") {
      return(NULL)
    }
    attempts <- diagnostics$attempts[[index]]
    if (is.null(attempts) || NROW(attempts) == 0L) return(NULL)
    data.frame(
      replication = diagnostics$replication[index],
      condition = diagnostics$condition[index],
      method = diagnostics$method[index],
      attempt = seq_len(NROW(attempts)),
      attempts,
      row.names = NULL,
      check.names = FALSE
    )
  })
  mc_combine_rows(rows)
}

mc_mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

mc_median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

mc_ipcr_it_nonconvergence_flags <- function(failures) {
  if (NROW(failures) == 0L) return(logical())
  flags <- rep(FALSE, NROW(failures))
  if ("ipcr_it_nonconvergence" %in% names(failures)) {
    explicit <- failures$ipcr_it_nonconvergence
    flags <- flags | (!is.na(explicit) & as.logical(explicit))
  }
  # This also recognizes checkpoints written after error classes began being
  # retained but before the explicit indicator was added.
  if ("error_classes" %in% names(failures)) {
    classes <- failures$error_classes
    flags <- flags | (!is.na(classes) & grepl(
      "ipcr_it_nonconvergence", classes, fixed = TRUE
    ))
  }
  flags
}

mc_summarize_diagnostics <- function(
    diagnostics, failures = NULL, settings = NULL,
    completed_replications = NULL) {
  if (is.null(failures)) failures <- mc_empty_data_frame()
  diagnostic_keys <- if (NROW(diagnostics) == 0L) {
    data.frame(condition = character(), method = character())
  } else {
    diagnostics[c("condition", "method")]
  }
  named_failures <- if (NROW(failures) == 0L ||
      !all(c("condition", "method") %in% names(failures))) {
    failures[FALSE, , drop = FALSE]
  } else {
    failures[
      !is.na(failures$condition) & nzchar(failures$condition) &
        !is.na(failures$method) & nzchar(failures$method),
      ,
      drop = FALSE
    ]
  }
  failure_keys <- if (NROW(named_failures) == 0L) {
    data.frame(condition = character(), method = character())
  } else {
    named_failures[c("condition", "method")]
  }
  planned_keys <- if (!is.null(settings) &&
      length(settings$conditions) > 0L && length(settings$methods) > 0L) {
    expand.grid(
      condition = as.character(settings$conditions),
      method = as.character(settings$methods),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(condition = character(), method = character())
  }
  group_keys <- unique(rbind(planned_keys, diagnostic_keys, failure_keys))
  if (NROW(group_keys) == 0L) return(mc_empty_data_frame())
  group_keys <- group_keys[
    order(group_keys$condition, group_keys$method),
    ,
    drop = FALSE
  ]

  mc_combine_rows(lapply(seq_len(NROW(group_keys)), function(index) {
    condition <- group_keys$condition[index]
    method <- group_keys$method[index]
    group <- diagnostics[
      diagnostics$condition == condition & diagnostics$method == method,
      ,
      drop = FALSE
    ]
    group_failures <- named_failures[
      named_failures$condition == condition &
        named_failures$method == method,
      ,
      drop = FALSE
    ]
    returned_replication_ids <- if (NROW(group) == 0L) integer() else {
      unique(as.integer(group$replication))
    }
    failure_replication_ids <- if (NROW(group_failures) == 0L) integer() else {
      unique(as.integer(group_failures$replication))
    }
    hard_failure_ids <- setdiff(
      failure_replication_ids, returned_replication_ids
    )
    nonconvergence_failure_ids <- if (NROW(group_failures) == 0L) {
      integer()
    } else {
      unique(as.integer(group_failures$replication[
        mc_ipcr_it_nonconvergence_flags(group_failures)
      ]))
    }
    attempted_replication_ids <- union(
      returned_replication_ids, failure_replication_ids
    )
    is_iterated <- identical(method, "iterated")
    converged_ids <- if (!is_iterated || NROW(group) == 0L) integer() else {
      unique(as.integer(group$replication[
        !is.na(group$converged) & group$converged
      ]))
    }
    nonconverged_ids <- if (!is_iterated || NROW(group) == 0L) {
      integer()
    } else {
      unique(as.integer(group$replication[
        !is.na(group$converged) & !group$converged
      ]))
    }
    unknown_status_ids <- if (!is_iterated || NROW(group) == 0L) {
      integer()
    } else {
      unique(as.integer(group$replication[is.na(group$converged)]))
    }
    attempted_count <- length(attempted_replication_ids)
    intended_count <- if (is.null(completed_replications)) {
      NA_integer_
    } else {
      as.integer(length(unique(completed_replications)))
    }
    convergence_rate_attempted <- if (is_iterated && attempted_count > 0L) {
      length(converged_ids) / attempted_count
    } else {
      NA_real_
    }
    convergence_rate_completed <- if (is_iterated &&
        !is.na(intended_count) && intended_count > 0L) {
      length(converged_ids) / intended_count
    } else {
      NA_real_
    }
    C_SE <- if ("C_SE" %in% names(group)) group$C_SE else numeric()
    C_l <- if ("C_l" %in% names(group)) group$C_l else numeric()
    data.frame(
      condition = condition,
      method = method,
      successful_replications = if (is_iterated) {
        length(converged_ids)
      } else {
        length(returned_replication_ids)
      },
      returned_replications = length(returned_replication_ids),
      converged_replications = if (is_iterated) {
        length(converged_ids)
      } else {
        NA_integer_
      },
      legacy_nonconverged_result_replications = if (is_iterated) {
        length(nonconverged_ids)
      } else {
        NA_integer_
      },
      unknown_convergence_status_replications = if (is_iterated) {
        length(unknown_status_ids)
      } else {
        NA_integer_
      },
      hard_analysis_failure_replications = length(hard_failure_ids),
      ipcr_it_nonconvergence_replications = if (is_iterated) {
        length(intersect(hard_failure_ids, nonconvergence_failure_ids))
      } else {
        NA_integer_
      },
      other_hard_analysis_failure_replications = if (is_iterated) {
        length(setdiff(hard_failure_ids, nonconvergence_failure_ids))
      } else {
        length(hard_failure_ids)
      },
      attempted_replications = attempted_count,
      completed_design_replications = intended_count,
      not_attempted_after_upstream_failure = if (is.na(intended_count)) {
        NA_integer_
      } else {
        max(0L, intended_count - attempted_count)
      },
      mean_pooled_fit_seconds = if (NROW(group) == 0L) {
        NA_real_
      } else {
        mc_mean_or_na(group$pooled_fit_seconds)
      },
      mean_elapsed_seconds = if (NROW(group) == 0L) {
        NA_real_
      } else {
        mc_mean_or_na(group$elapsed_seconds)
      },
      convergence_rate = convergence_rate_attempted,
      convergence_rate_among_attempted = convergence_rate_attempted,
      convergence_rate_among_completed_replications =
        convergence_rate_completed,
      mean_iterations = if (is_iterated && NROW(group) > 0L) {
        mc_mean_or_na(group$iterations)
      } else {
        NA_real_
      },
      mean_total_iterations = if (is_iterated && NROW(group) > 0L) {
        mc_mean_or_na(group$total_iterations)
      } else {
        NA_real_
      },
      mean_learning_rate_attempts = if (is_iterated && NROW(group) > 0L) {
        mc_mean_or_na(group$learning_rate_attempts)
      } else {
        NA_real_
      },
      mean_selected_learning_rate = if (is_iterated && NROW(group) > 0L) {
        mc_mean_or_na(group$learning_rate)
      } else {
        NA_real_
      },
      median_C_SE = if (is_iterated) mc_median_or_na(C_SE) else NA_real_,
      median_C_l = if (is_iterated) mc_median_or_na(C_l) else NA_real_,
      replications_with_fit_warnings = if (NROW(group) == 0L) {
        0L
      } else {
        sum(nzchar(group$pooled_fit_warnings))
      },
      replications_with_warnings = if (NROW(group) == 0L) {
        0L
      } else {
        sum(nzchar(group$warnings))
      },
      row.names = NULL
    )
  }))
}

mc_pillai_test <- function(result, term) {
  if (is.null(result$mancova)) {
    return(c(statistic = NA_real_, df1 = NA_real_, df2 = NA_real_,
             p_value = NA_real_))
  }
  mancova_summary <- summary(
    result$mancova,
    test.statistic = "Pillai",
    multivariate = TRUE,
    univariate = FALSE
  )
  test <- mancova_summary$multivariate.tests[[term]]
  if (is.null(test) || isTRUE(test$singular)) {
    return(c(statistic = NA_real_, df1 = NA_real_, df2 = NA_real_,
             p_value = NA_real_))
  }

  eigenvalues <- Re(eigen(
    qr.coef(qr(test$SSPE), test$SSPH),
    symmetric = FALSE,
    only.values = TRUE
  )$values)
  statistic <- sum(eigenvalues / (1 + eigenvalues))
  p <- length(eigenvalues)
  q <- unname(test$df)
  s <- min(p, q)
  n <- 0.5 * (test$df.residual - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  temporary_1 <- 2 * m + s + 1
  temporary_2 <- 2 * n + s + 1
  approximate_f <- (temporary_2 / temporary_1 * statistic) /
    (s - statistic)
  df1 <- s * temporary_1
  df2 <- s * temporary_2
  p_value <- stats::pf(
    approximate_f,
    df1 = df1,
    df2 = df2,
    lower.tail = FALSE
  )
  c(statistic = statistic, df1 = df1, df2 = df2, p_value = p_value)
}

mc_summarize_tests <- function(tests) {
  if (NROW(tests) == 0L) return(mc_empty_data_frame())
  groups <- split(
    tests,
    interaction(tests$condition, tests$method, tests$test, drop = TRUE)
  )
  mc_combine_rows(lapply(groups, function(group) {
    valid <- is.finite(group$p_value)
    group <- group[valid, , drop = FALSE]
    if (NROW(group) == 0L) return(NULL)
    data.frame(
      condition = group$condition[1L],
      method = group$method[1L],
      test = group$test[1L],
      rejection_rate_05 = mean(group$p_value < 0.05),
      successful_replications = NROW(group),
      row.names = NULL
    )
  }))
}

mc_checkpoint_payload <- function(results, settings, elapsed = NULL) {
  estimates <- mc_combine_rows(lapply(results, function(x) x$estimates))
  diagnostics <- mc_combine_rows(lapply(results, function(x) x$diagnostics))
  tests <- mc_combine_rows(lapply(results, function(x) x$tests))
  failures <- mc_combine_rows(lapply(results, function(x) x$failures))
  # Preserve all returned rows for diagnosis, but never use a legacy
  # nonconverged iterated-IPCR result in an inferential simulation summary.
  converged_estimates <- mc_keep_converged_iterated(estimates, diagnostics)
  converged_tests <- mc_keep_converged_iterated(tests, diagnostics)
  completed_replications <- which(!vapply(results, is.null, logical(1)))
  list(
    settings = settings,
    estimates = estimates,
    diagnostics = diagnostics,
    attempt_diagnostics = mc_combine_rows(list(
      mc_collect_attempt_diagnostics(diagnostics),
      mc_collect_attempt_diagnostics(failures)
    )),
    tests = tests,
    failures = failures,
    summary = mc_summarize_estimates(converged_estimates),
    diagnostic_summary = mc_summarize_diagnostics(
      diagnostics, failures, settings, completed_replications
    ),
    test_summary = mc_summarize_tests(converged_tests),
    elapsed = elapsed,
    saved_at = Sys.time()
  )
}

mc_write_checkpoint <- function(results, settings, path, elapsed = NULL) {
  payload <- mc_checkpoint_payload(results, settings, elapsed)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(payload, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop(
      "Unable to atomically replace the Monte Carlo checkpoint at '",
      path, "'.",
      call. = FALSE
    )
  }
  payload
}

mc_elapsed_since <- function(started_at, previous_seconds = 0) {
  elapsed_seconds <- previous_seconds + as.numeric(
    difftime(Sys.time(), started_at, units = "secs")
  )
  as.difftime(elapsed_seconds, units = "secs")
}

mc_resume_checkpoint <- function(
    path, settings, n_replications,
    ignored_settings = c("n_replications", "checkpoint_every")) {
  empty_resume <- function() {
    list(
      results = vector("list", n_replications),
      completed_replications = 0L,
      next_replication = 1L,
      counts = c(failures = 0L, nonconvergence_errors = 0L),
      elapsed_seconds = 0
    )
  }
  if (!file.exists(path)) return(empty_resume())

  checkpoint <- tryCatch(
    readRDS(path),
    error = function(error) {
      stop(
        "Unable to read the existing Monte Carlo checkpoint at '", path,
        "': ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  required_components <- c(
    "settings", "estimates", "diagnostics", "tests", "failures"
  )
  missing_components <- setdiff(required_components, names(checkpoint))
  if (length(missing_components) > 0L) {
    stop(
      "The existing Monte Carlo checkpoint is missing: ",
      paste(missing_components, collapse = ", "), ".",
      call. = FALSE
    )
  }

  settings_to_check <- setdiff(names(settings), ignored_settings)
  incompatible_settings <- settings_to_check[vapply(
    settings_to_check,
    function(setting) {
      !setting %in% names(checkpoint$settings) ||
        !isTRUE(all.equal(
          checkpoint$settings[[setting]], settings[[setting]],
          check.attributes = FALSE
        ))
    },
    FUN.VALUE = logical(1)
  )]
  if (length(incompatible_settings) > 0L) {
    stop(
      "The existing Monte Carlo checkpoint is incompatible with the current ",
      "settings: ", paste(incompatible_settings, collapse = ", "),
      ". Use a different results file or remove the incompatible checkpoint.",
      call. = FALSE
    )
  }

  result_components <- checkpoint[c(
    "estimates", "diagnostics", "tests", "failures"
  )]
  completed_replications <- sort(unique(unlist(lapply(
    result_components,
    function(component) {
      if (NROW(component) == 0L ||
          !"replication" %in% names(component)) {
        return(integer())
      }
      as.integer(component$replication)
    }
  ))))
  if (length(completed_replications) == 0L) return(empty_resume())
  if (anyNA(completed_replications) || any(completed_replications < 1L)) {
    stop(
      "The existing Monte Carlo checkpoint has invalid replication indices.",
      call. = FALSE
    )
  }

  last_replication <- max(completed_replications)
  if (!identical(completed_replications, seq_len(last_replication))) {
    stop(
      "The existing Monte Carlo checkpoint does not contain a contiguous ",
      "sequence of completed replications.",
      call. = FALSE
    )
  }
  if (last_replication > n_replications) {
    stop(
      "The checkpoint already contains ", last_replication,
      " replications, more than the requested ", n_replications, ".",
      call. = FALSE
    )
  }

  replication_results <- vector("list", n_replications)
  for (replication in completed_replications) {
    replication_results[[replication]] <- lapply(
      result_components,
      function(component) {
        if (NROW(component) == 0L ||
            !"replication" %in% names(component)) {
          return(mc_empty_data_frame())
        }
        component[
          component$replication == replication,
          ,
          drop = FALSE
        ]
      }
    )
  }

  cumulative_counts <- Reduce(
    `+`,
    lapply(
      replication_results[completed_replications],
      mc_replication_counts
    ),
    init = c(failures = 0L, nonconvergence_errors = 0L)
  )
  elapsed_seconds <- if (is.null(checkpoint$elapsed)) {
    0
  } else if (inherits(checkpoint$elapsed, "difftime")) {
    as.numeric(checkpoint$elapsed, units = "secs")
  } else {
    as.numeric(checkpoint$elapsed)
  }
  if (length(elapsed_seconds) != 1L || !is.finite(elapsed_seconds) ||
      elapsed_seconds < 0) {
    stop(
      "The existing Monte Carlo checkpoint has an invalid elapsed time.",
      call. = FALSE
    )
  }

  list(
    results = replication_results,
    completed_replications = last_replication,
    next_replication = last_replication + 1L,
    counts = cumulative_counts,
    elapsed_seconds = elapsed_seconds
  )
}

mc_failure <- function(
    replication, stage, message, condition = NA_character_,
    method = NA_character_, error = NULL) {
  is_iterated_nonconvergence <- !is.null(error) &&
    inherits(error, "ipcr_it_nonconvergence")
  attempts <- if (is_iterated_nonconvergence &&
      !is.null(error$attempts)) {
    error$attempts
  } else {
    data.frame()
  }
  iteration_history <- if (is_iterated_nonconvergence) {
    error$iteration_history
  } else {
    NULL
  }
  list(
    estimates = mc_empty_data_frame(),
    diagnostics = mc_empty_data_frame(),
    tests = mc_empty_data_frame(),
    failures = data.frame(
      replication = replication,
      condition = condition,
      method = method,
      stage = stage,
      message = message,
      ipcr_it_nonconvergence = is_iterated_nonconvergence,
      error_classes = if (is.null(error)) {
        NA_character_
      } else {
        paste(class(error), collapse = " | ")
      },
      attempts = I(list(attempts)),
      iteration_history = I(list(iteration_history)),
      row.names = NULL
    )
  )
}
