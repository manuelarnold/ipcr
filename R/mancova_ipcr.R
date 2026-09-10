fit_mancova_ipcr <- function(mlm, q) {
  mancova <- tryCatch(
    car::Manova(
      mod = mlm,
      type = "II",
      test.statistic = "Pillai"
    ),
    error = function(e) {
      warning(
        "Type-II MANCOVA is unavailable: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  if (is.null(mancova)) {
    return(list(model = NULL, summary = NULL))
  }

  residual_ssp <- mancova$SSPE
  residual_rank <- if (all(is.finite(residual_ssp))) {
    qr(residual_ssp)$rank
  } else {
    0L
  }
  if (residual_rank < q) {
    warning(
      "Type-II MANCOVA is unavailable because the residual IPC matrix has ",
      "rank ", residual_rank, " but contains ", q, " responses. The ",
      "individual IPC regressions and their F-tests remain available.",
      call. = FALSE
    )
    return(list(model = mancova, summary = NULL))
  }

  mancova_summary <- tryCatch(
    summary(
      mancova,
      test.statistic = "Pillai",
      multivariate = TRUE,
      univariate = FALSE
    ),
    error = function(e) {
      warning(
        "Type-II MANCOVA summary is unavailable: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  # car computes the multivariate statistics in its print method. Validate the
  # stored summary here so that a singular or otherwise invalid result is
  # reported by ipcr(), rather than later by summary().
  if (!is.null(mancova_summary)) {
    valid_summary <- tryCatch(
      {
        utils::capture.output(
          print(mancova_summary, SSP = FALSE, SSPE = FALSE)
        )
        TRUE
      },
      error = function(e) {
        warning(
          "Type-II MANCOVA summary is unavailable: ",
          conditionMessage(e),
          call. = FALSE
        )
        FALSE
      }
    )
    if (!valid_summary) {
      mancova_summary <- NULL
    }
  }

  list(model = mancova, summary = mancova_summary)
}
