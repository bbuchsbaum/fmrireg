make_subset_contrast_fixture <- function(seed = 1L) {
  set.seed(seed)
  n_trials <- 20L
  trial_on <- seq(5, by = 20, length.out = n_trials)

  events <- data.frame(
    run = 1L,
    onset = as.vector(rbind(trial_on, trial_on + 4, trial_on + 10)),
    duration = rep(c(4, 6, 2), times = n_trials),
    phase = factor(
      rep(c("a", "b", "c"), times = n_trials),
      levels = c("a", "b", "c")
    ),
    modulator = rep(scale(rnorm(n_trials), scale = FALSE)[, 1], each = 3)
  )
  events$phase_pm <- factor(as.character(events$phase), levels = c("a", "b"))

  list(
    events = events,
    sframe = sampling_frame(blocklens = 250, TR = 2),
    response = matrix(rnorm(250 * 5), 250, 5)
  )
}

make_subset_contrast_model <- function(fixture, workaround = FALSE) {
  events <- fixture$events
  sframe <- fixture$sframe

  if (workaround) {
    pm_contrasts <- contrast_set(
      column_contrast("^phase_pm[.]a_", name = "pm_a"),
      column_contrast("^phase_pm[.]b_", name = "pm_b")
    )
    event <- event_model(
      onset ~ hrf(phase, basis = "spmg1") +
        hrf(
          phase_pm, modulator,
          basis = "spmg1", prefix = "pm",
          subset = !is.na(phase_pm), contrasts = pm_contrasts
        ),
      data = events, block = ~run, sampling_frame = sframe,
      durations = events$duration
    )
  } else {
    pm_contrasts <- contrast_set(
      column_contrast("^phase[.]a_", name = "pm_a"),
      column_contrast("^phase[.]b_", name = "pm_b")
    )
    event <- event_model(
      onset ~ hrf(phase, basis = "spmg1") +
        hrf(
          phase, modulator,
          basis = "spmg1", prefix = "pm",
          subset = phase != "c", contrasts = pm_contrasts
        ),
      data = events, block = ~run, sampling_frame = sframe,
      durations = events$duration
    )
  }

  dataset <- matrix_dataset(
    fixture$response,
    TR = 2, run_length = 250, event_table = events
  )
  fmri_model(
    event,
    baseline_model(basis = "bs", degree = 3, sframe = sframe),
    dataset
  )
}

fit_subset_contrast_model <- function(model, fast) {
  warnings <- character(0)
  fit <- withCallingHandlers(
    fmri_lm(
      model,
      strategy = "chunkwise", nchunks = 1,
      use_fast_path = fast, progress = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warnings)
}

test_that("subsetted hrf t-contrasts use reconciled design-wide weights", {
  fixture <- make_subset_contrast_fixture()
  subset_model <- make_subset_contrast_model(fixture)
  workaround_model <- suppressWarnings(
    make_subset_contrast_model(fixture, workaround = TRUE)
  )

  # The workaround changes only the factor used to name the two retained
  # modulator columns; the numeric design is an oracle for the subset model.
  expect_equal(
    unname(design_matrix(subset_model)),
    unname(design_matrix(workaround_model)),
    tolerance = 1e-12
  )

  reference <- fit_subset_contrast_model(workaround_model, fast = TRUE)$fit
  reference_stats <- as.matrix(stats(reference, type = "contrasts"))

  for (fast in c(TRUE, FALSE)) {
    observed <- fit_subset_contrast_model(subset_model, fast = fast)

    expect_setequal(observed$fit$result$contrasts$name, c("pm_a", "pm_b"))
    expect_false(any(grepl("Failed to compute t-contrast", observed$warnings)))
    expect_equal(
      unname(as.matrix(stats(observed$fit, type = "contrasts"))),
      unname(reference_stats),
      tolerance = 1e-8
    )
  }
})

test_that("invalid requested contrasts fail instead of being dropped", {
  bad <- matrix(c(1, 0, 0), ncol = 1)
  attr(bad, "colind") <- c(1L, 2L)

  expect_error(
    process_t_contrasts(
      B = matrix(0, 2, 1), sigma2 = 1, XtXinv = diag(2),
      conlist = list(requested = bad), df = 8
    ),
    "Failed to compute t-contrast 'requested'.*Dimension mismatch"
  )

  bad_f <- matrix(1, nrow = 2, ncol = 2)
  attr(bad_f, "colind") <- c(1L, 2L, 3L)
  expect_error(
    process_f_contrasts(
      B = matrix(0, 3, 1), sigma2 = 1, XtXinv = diag(3),
      fconlist = list(requested_F = bad_f), df = 8
    ),
    "Failed to compute F-contrast 'requested_F'.*do not match"
  )
})
