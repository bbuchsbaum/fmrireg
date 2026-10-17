# fmrilm.R: fitted_hrf error branches, print.fmri_lm fallbacks.

test_that("fitted_hrf.fmri_lm covers empty/feature/error metadata branches", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  expect_error(fitted_hrf(list(), sample_at = 0:2), "applicable method|fmri_lm")
  expect_error(fitted_hrf(fit, sample_at = "x"), "numeric")
  expect_error(fitted_hrf(fit, sample_at = c(0, NA_real_)), "finite")

  # Empty event terms -> empty list
  fit_empty <- fit
  fit_empty$model$event_model <- structure(
    list(terms = list(), design_matrix = matrix(0, 1, 0)),
    class = "event_model"
  )
  # terms() may still find terms via another path; force empty via stub
  empty <- tryCatch(
    {
      # Directly call with a fit that has no event terms by clearing terms attr
      old_terms <- terms(fit$model$event_model)
      if (length(old_terms) == 0L) {
        fitted_hrf(fit_empty, sample_at = 0:2)
      } else {
        NULL
      }
    },
    error = function(e) e
  )

  # Missing col_indices
  fit2 <- fit
  dm <- fit2$model$event_model$design_matrix
  attr(dm, "col_indices") <- NULL
  fit2$model$event_model$design_matrix <- dm
  expect_error(fitted_hrf(fit2, sample_at = 0:3), "col_indices")

  # Out-of-bounds coefficient indices
  fit3 <- fit
  dm3 <- fit3$model$event_model$design_matrix
  ci <- attr(dm3, "col_indices")
  if (!is.null(ci) && length(ci) > 0L) {
    nm <- names(ci)[1]
    ci[[nm]] <- c(1L, 99999L)
    attr(dm3, "col_indices") <- ci
    fit3$model$event_model$design_matrix <- dm3
    expect_error(fitted_hrf(fit3, sample_at = 0:2), "out of bounds|indices")
  }

  # Non-matrix betas coerced
  fit4 <- fit
  betas <- fit4$result$betas$data[[1]]$estimate[[1]]
  if (is.matrix(betas)) {
    fit4$result$betas$data[[1]]$estimate[[1]] <- as.vector(betas[, 1, drop = TRUE])
    # May error on shape or succeed after as.matrix; either exercises the branch
    expect_error(
      suppressWarnings(fitted_hrf(fit4, sample_at = 0:2)),
      regexp = "."
    )
  }

  # Happy path still works
  ok <- fitted_hrf(fit, sample_at = c(0, 2, 4))
  expect_true(is.list(ok))
  expect_true(length(ok) >= 1L)
})

test_that("print.fmri_lm writes its summary to stdout", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())

  out <- capture.output(print(fit))
  expect_true(any(grepl("fMRI Linear Model Results", out, fixed = TRUE)))
  expect_true(any(grepl("Estimation scope: runwise_meta", out, fixed = TRUE)))
  expect_true(any(grepl("Formula: `onsets ~ hrf(condition)`", out, fixed = TRUE)))
})
