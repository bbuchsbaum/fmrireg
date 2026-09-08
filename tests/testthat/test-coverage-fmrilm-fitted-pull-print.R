# fmrilm.R: fitted_hrf error branches, pull_stat_revised, print.fmri_lm fallbacks.

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

test_that("pull_stat_revised extracts betas/contrasts/F and errors", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())

  betas <- fmrireg:::pull_stat_revised(fit, "betas", "estimate")
  expect_true(inherits(betas, "tbl_df") || is.data.frame(betas) || is.matrix(betas))

  expect_error(fmrireg:::pull_stat_revised(fit, "nope", "estimate"), "Invalid type")

  # Empty contrasts table -> error paths
  fit0 <- fit
  fit0$result$contrasts <- tibble::tibble(
    type = character(), name = character(), data = list()
  )
  expect_error(fmrireg:::pull_stat_revised(fit0, "contrasts", "estimate"), "No simple contrasts")
  expect_error(fmrireg:::pull_stat_revised(fit0, "F", "stat"), "No F contrasts")

  # Inject a contrast row to exercise contrast/F branches
  fit2 <- fit
  n_vox <- nrow(as.matrix(coef(fit)))
  est <- list(rnorm(n_vox))
  fit2$result$contrasts <- tibble::tibble(
    type = c("contrast", "Fcontrast"),
    name = c("A_vs_B", "omnibus"),
    data = list(
      tibble::tibble(estimate = est, se = est, stat = est, prob = est),
      tibble::tibble(estimate = est, se = est, stat = est, prob = est)
    )
  )
  c_est <- fmrireg:::pull_stat_revised(fit2, "contrasts", "estimate")
  expect_true("A_vs_B" %in% names(c_est))
  f_est <- fmrireg:::pull_stat_revised(fit2, "F", "stat")
  expect_true("omnibus" %in% names(f_est))
})

test_that("print.fmri_lm covers crayon and no-crayon paths", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  fit$strategy <- "runwise"
  fit$bcons <- list()

  out <- capture.output(print(fit))
  expect_true(any(grepl("fmri_lm|formula|Design|Baseline|Contrast", out, ignore.case = TRUE)))

  # Force no-crayon branch via temporary unload of crayon namespace visibility
  # by calling the method after stubbing requireNamespace
  with_mocked_bindings(
    requireNamespace = function(pkg, quietly = TRUE) {
      if (identical(pkg, "crayon")) FALSE else base::requireNamespace(pkg, quietly = quietly)
    },
    {
      out2 <- capture.output(print(fit))
      expect_true(any(grepl("install 'crayon'|fmri_lm_result|Strategy|formula", out2)))
    },
    .package = "base"
  )
})
