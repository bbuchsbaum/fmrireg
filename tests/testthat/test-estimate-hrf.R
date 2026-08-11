.make_estimate_hrf_fixture <- function(
    n = 220L,
    n_voxels = 4L,
    noise_sd = 0.025,
    include_fixed = TRUE,
    seed = 180L) {
  set.seed(seed)
  onsets <- seq(8, n - 30, by = 10)
  events <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = length(onsets))),
    nuisance = sin(seq_along(onsets) * 1.7),
    run = 1L
  )
  sampling <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  event_model <- fmridesign::event_model(
    onset ~ hrf(condition),
    data = events,
    block = ~run,
    sampling_frame = sampling
  )
  X_event <- as.matrix(fmridesign::design_matrix(event_model))

  amplitude_a <- rep(c(1.0, 0.7, -0.6, 1.4), length.out = n_voxels)
  amplitudes <- rbind(A = amplitude_a, B = 2 * amplitude_a)
  Y <- X_event %*% amplitudes

  if (include_fixed) {
    fixed_model <- fmridesign::event_model(
      onset ~ hrf(nuisance),
      data = events,
      block = ~run,
      sampling_frame = sampling
    )
    X_fixed <- as.matrix(fmridesign::design_matrix(fixed_model))
    Y <- Y + X_fixed %*% matrix(seq(0.8, 1.1, length.out = n_voxels), nrow = 1L)
  }

  Y <- sweep(Y, 2L, seq(-0.3, 0.3, length.out = n_voxels), "+")
  Y <- Y + matrix(rnorm(n * n_voxels, sd = noise_sd), nrow = n)
  colnames(Y) <- paste0("voxel_", seq_len(n_voxels))

  list(
    dataset = fmridataset::matrix_dataset(
      datamat = Y,
      TR = 1,
      run_length = n,
      event_table = events
    ),
    events = events,
    Y = Y,
    amplitudes = amplitudes
  )
}

.make_delayed_hrf_fixture <- function(seed = 731L) {
  set.seed(seed)
  n <- 300L
  onsets <- 8 + c(0, cumsum(sample(c(6, 7, 9, 10), 31L, replace = TRUE)))
  events <- data.frame(
    onset = onsets,
    condition = factor(sample(rep(c("A", "B"), 16L))),
    run = 1L
  )
  sampling <- fmrihrf::sampling_frame(blocklens = n, TR = 1)
  acquisition_time <- fmrihrf::samples(sampling, global = TRUE)
  truth_a <- function(time) as.numeric(fmrihrf::HRF_SPMG1(time))
  truth_b <- function(time) {
    ifelse(
      time >= 1.5,
      1.25 * as.numeric(fmrihrf::HRF_SPMG1(time - 1.5)),
      0
    )
  }
  signal_a <- fmrihrf::evaluate(
    fmrihrf::regressor(
      events$onset[events$condition == "A"], hrf = truth_a, span = 24
    ),
    acquisition_time
  )
  signal_b <- fmrihrf::evaluate(
    fmrihrf::regressor(
      events$onset[events$condition == "B"], hrf = truth_b, span = 24
    ),
    acquisition_time
  )
  amplitudes <- rbind(A = c(1, 0.65), B = c(0.8, 1.1))
  drift <- outer(seq(-0.15, 0.15, length.out = n), rep(1, 2L))
  response <- signal_a %o% amplitudes["A", ] +
    signal_b %o% amplitudes["B", ] + drift +
    matrix(rnorm(n * 2L, sd = 0.025), nrow = n)
  colnames(response) <- c("visual_ROI", "motor_ROI")
  dataset <- fmridataset::matrix_dataset(
    response, TR = 1, run_length = n, event_table = events
  )

  list(
    dataset = dataset,
    baseline = fmridesign::baseline_model(
      "poly", degree = 1, sframe = dataset$sampling_frame
    ),
    amplitudes = amplitudes,
    truth = list(A = truth_a, B = truth_b)
  )
}

test_that("estimate_hrf recovers labeled condition curves with fixed effects", {
  fixture <- .make_estimate_hrf_fixture()
  sample_at <- 0:20

  fit <- estimate_hrf(
    onset ~ hrf(condition),
    fixed = onset ~ hrf(nuisance),
    block = ~run,
    dataset = fixture$dataset,
    rsam = sample_at,
    k = 8L,
    lambda = "gcv",
    progress = FALSE
  )

  expect_s3_class(fit, "fmri_hrf_estimate")
  expect_identical(dim(fit$estimate), c(length(sample_at), 2L, 4L))
  expect_identical(dim(fit$std.error), dim(fit$estimate))
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$std.error)))
  expect_true(all(fit$std.error >= 0))
  expect_true(fit$lambda %in% fit$gcv$lambda)
  expect_true(all(c("curve", "term", "condition") %in% names(fit$curve_info)))

  truth_shape <- as.numeric(fmrihrf::HRF_SPMG1(sample_at))
  for (curve_index in seq_len(nrow(fit$curve_info))) {
    condition <- sub("^.*\\.", "", fit$curve_info$condition[[curve_index]])
    truth_row <- match(condition, rownames(fixture$amplitudes))
    expect_false(is.na(truth_row))
    for (voxel in seq_len(ncol(fixture$amplitudes))) {
      truth <- truth_shape * fixture$amplitudes[truth_row, voxel]
      observed <- fit$estimate[, curve_index, voxel]
      expect_gt(stats::cor(observed, truth), 0.9)
      expect_lt(sqrt(mean((observed - truth)^2)) / max(abs(truth)), 0.35)
    }
  }

  expect_true(all(abs(fit$estimate[1L, , ]) < 1e-12))
  expect_true(all(abs(fit$estimate[length(sample_at), , ]) < 1e-12))
})

test_that("estimate_hrf distinguishes canonical and delayed response shapes", {
  fixture <- .make_delayed_hrf_fixture()
  fit <- estimate_hrf(
    onset ~ hrf(condition),
    block = ~run,
    dataset = fixture$dataset,
    basemod = fixture$baseline,
    rsam = seq(0, 24, by = 0.5),
    k = 8L,
    lambda = "gcv"
  )
  peak_times <- matrix(
    NA_real_, nrow = 2L, ncol = length(fit$voxels),
    dimnames = list(c("A", "B"), fit$voxels)
  )

  for (curve_index in seq_along(fit$curves)) {
    condition <- sub("^.*[.]", "", fit$curve_info$condition[[curve_index]])
    truth_shape <- fixture$truth[[condition]](fit$time)
    for (voxel in seq_along(fit$voxels)) {
      truth <- truth_shape * fixture$amplitudes[condition, voxel]
      observed <- fit$estimate[, curve_index, voxel]
      peak_times[condition, voxel] <- fit$time[[which.max(observed)]]

      expect_gt(stats::cor(observed, truth), 0.97)
      expect_lt(
        sqrt(mean((observed - truth)^2)) / max(abs(truth)),
        0.08
      )
    }
  }

  expect_true(all(peak_times["B", ] - peak_times["A", ] >= 1.5))
  expect_true(is.finite(fit$lambda) && fit$lambda >= 0)
})

test_that("estimate_hrf is linear in the response for fixed smoothing", {
  fixture <- .make_estimate_hrf_fixture(n_voxels = 2L, include_fixed = FALSE)
  scaled_dataset <- fmridataset::matrix_dataset(
    datamat = 3.5 * fixture$Y,
    TR = 1,
    run_length = nrow(fixture$Y),
    event_table = fixture$events
  )

  fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
    rsam = 0:20, k = 8L, lambda = 1
  )
  scaled <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = scaled_dataset,
    rsam = 0:20, k = 8L, lambda = 1
  )

  expect_equal(scaled$estimate, 3.5 * fit$estimate, tolerance = 1e-9)
  expect_equal(scaled$std.error, 3.5 * fit$std.error, tolerance = 1e-9)
  expect_equal(scaled$lambda, fit$lambda)
})

test_that("estimate_hrf methods expose curves without losing labels", {
  fixture <- .make_estimate_hrf_fixture(n_voxels = 2L, include_fixed = FALSE)
  fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
    rsam = 0:12, k = 6L, lambda = 0.5
  )

  tbl <- tidy(fit)
  expect_s3_class(tbl, "tbl_df")
  expect_identical(
    names(tbl),
    c("time", "curve", "term", "condition", "voxel", "estimate", "std.error", "lower", "upper")
  )
  expect_equal(nrow(tbl), length(fit$time) * length(fit$curves) * length(fit$voxels))

  first_curve <- fit$curves[[1L]]
  expect_equal(as.matrix(fit, curve = first_curve), fit$estimate[, first_curve, ])
  prediction <- predict(fit, newdata = seq(0, fit$span, length.out = 17L), se.fit = TRUE)
  expect_identical(dim(prediction$fit), c(17L, 2L, 2L))
  expect_identical(dim(prediction$se.fit), dim(prediction$fit))
  expect_true(all(is.finite(prediction$fit)))
})

test_that("estimate_hrf supports tent and explicit legacy compatibility paths", {
  fixture <- .make_estimate_hrf_fixture(n_voxels = 1L, include_fixed = FALSE)

  tent <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
    rsam = 0:12, basis = "tent", k = 5L, lambda = 1
  )
  expect_identical(tent$basis, "tent")
  expect_true(all(is.finite(tent$estimate)))
  expect_true(all(abs(tent$estimate[c(1L, length(tent$time)), , ]) < 1e-12))

  expect_warning(
    legacy_bs <- estimate_hrf(
      onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
      bs = "tp", rsam = 0:12, k = 6L, lambda = 1
    ),
    "bs is deprecated"
  )
  expect_identical(legacy_bs$basis, "bspline")

  expect_warning(
    legacy_fx <- estimate_hrf(
      onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
      rsam = 0:12, k = 6L, fx = TRUE
    ),
    "fx is deprecated"
  )
  expect_identical(legacy_fx$lambda, 0)
})

test_that("estimate_hrf validates the estimand and edge cases", {
  fixture <- .make_estimate_hrf_fixture(n_voxels = 1L, include_fixed = FALSE)

  expect_error(
    estimate_hrf(onset ~ trialwise(), block = ~run, dataset = fixture$dataset),
    "condition-level hrf"
  )
  expect_error(
    estimate_hrf(onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
                 rsam = c(0, 2, 1)),
    "strictly increasing"
  )
  expect_error(
    estimate_hrf(onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
                 rsam = 1:10),
    "start at zero"
  )
  expect_error(
    estimate_hrf(onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
                 k = 3L),
    "at least 4"
  )
  expect_error(
    predict(
      estimate_hrf(onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
                   rsam = 0:12, k = 6L, lambda = 1),
      newdata = c(-1, 2)
    ),
    "within"
  )
})

test_that("estimate_hrf uses one shared multiresponse solve", {
  fixture <- .make_estimate_hrf_fixture(
    n = 180L,
    n_voxels = 128L,
    noise_sd = 0.05,
    include_fixed = FALSE
  )

  elapsed <- system.time({
    fit <- estimate_hrf(
      onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
      rsam = 0:16, k = 6L, lambda = 1
    )
  })[["elapsed"]]

  expect_identical(fit$diagnostics$algorithm, "shared penalized multiresponse solve")
  expect_lt(elapsed, 10)
})

test_that("unpenalized estimate_hrf matches a direct full-model QR oracle", {
  fixture <- .make_estimate_hrf_fixture(
    n_voxels = 2L,
    noise_sd = 0.04,
    include_fixed = TRUE
  )
  fit <- estimate_hrf(
    onset ~ hrf(condition),
    fixed = onset ~ hrf(nuisance),
    block = ~run,
    dataset = fixture$dataset,
    rsam = 0:20,
    k = 8L,
    lambda = 0
  )

  full_design <- cbind(fit$nuisance_design, fit$event_design)
  full_qr <- qr(full_design, LAPACK = FALSE)
  direct_coefficients <- qr.coef(full_qr, fixture$Y)
  event_rows <- ncol(fit$nuisance_design) + seq_len(ncol(fit$event_design))
  direct_residuals <- fixture$Y - full_design %*% direct_coefficients
  direct_sigma <- sqrt(
    colSums(direct_residuals^2) / (nrow(full_design) - full_qr$rank)
  )
  direct_covariance <- chol2inv(chol(crossprod(full_design)))[
    event_rows, event_rows, drop = FALSE
  ]

  expect_equal(
    unname(fit$coefficients),
    unname(direct_coefficients[event_rows, , drop = FALSE]),
    tolerance = 1e-9
  )
  expect_equal(unname(fit$sigma), unname(direct_sigma), tolerance = 1e-9)
  expect_equal(
    unname(fit$covariance_unscaled),
    unname(direct_covariance),
    tolerance = 1e-8
  )
})

test_that("GCV and fitted curves respect scale and ordering invariants", {
  fixture <- .make_estimate_hrf_fixture(
    n_voxels = 3L,
    noise_sd = 0.04,
    include_fixed = FALSE
  )
  fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
    rsam = 0:20, k = 8L, lambda = "gcv"
  )

  scales <- c(0.1, 4, 25)
  scaled_data <- sweep(fixture$Y, 2L, scales, "*")
  scaled_dataset <- fmridataset::matrix_dataset(
    scaled_data,
    TR = 1,
    run_length = nrow(scaled_data),
    event_table = fixture$events
  )
  scaled_fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = scaled_dataset,
    rsam = 0:20, k = 8L, lambda = "gcv"
  )

  expect_equal(scaled_fit$gcv$score, fit$gcv$score, tolerance = 1e-10)
  expect_identical(scaled_fit$lambda, fit$lambda)
  expect_equal(
    sweep(scaled_fit$estimate, 3L, scales, "/"),
    fit$estimate,
    tolerance = 1e-9
  )

  permutation <- c(3L, 1L, 2L)
  permuted_dataset <- fmridataset::matrix_dataset(
    fixture$Y[, permutation, drop = FALSE],
    TR = 1,
    run_length = nrow(fixture$Y),
    event_table = fixture$events[sample.int(nrow(fixture$events)), , drop = FALSE]
  )
  permuted_fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = permuted_dataset,
    rsam = 0:20, k = 8L, lambda = "gcv"
  )
  original_order <- match(fit$voxels, permuted_fit$voxels)

  expect_identical(permuted_fit$lambda, fit$lambda)
  expect_equal(
    permuted_fit$estimate[, , original_order, drop = FALSE],
    fit$estimate,
    tolerance = 1e-9
  )
})

test_that("estimate_hrf handles runwise baselines and condition labels", {
  set.seed(812)
  run_length <- c(120L, 120L)
  events <- do.call(rbind, lapply(seq_along(run_length), function(run) {
    onsets <- 6 + c(0, cumsum(sample(c(6, 7, 9), 11L, replace = TRUE)))
    data.frame(
      onset = onsets,
      condition = factor(sample(rep(c("A", "B"), 6L))),
      run = run
    )
  }))
  sampling <- fmrihrf::sampling_frame(blocklens = run_length, TR = 1)
  event_model <- fmridesign::event_model(
    onset ~ hrf(condition),
    data = events,
    block = ~run,
    sampling_frame = sampling
  )
  X_event <- as.matrix(fmridesign::design_matrix(event_model))
  amplitudes <- rbind(A = c(1, 0.7), B = c(1.6, -0.8))
  Y <- X_event %*% amplitudes +
    rep(c(-0.4, 0.6), times = run_length) +
    matrix(rnorm(sum(run_length) * 2L, sd = 0.025), ncol = 2L)
  colnames(Y) <- c("left_ROI", "right_ROI")
  dataset <- fmridataset::matrix_dataset(
    Y, TR = 1, run_length = run_length, event_table = events
  )

  fit <- estimate_hrf(
    onset ~ hrf(condition), block = ~run, dataset = dataset,
    rsam = 0:20, k = 8L, lambda = "gcv"
  )
  truth_shape <- as.numeric(fmrihrf::HRF_SPMG1(fit$time))

  expect_identical(fit$diagnostics$nuisance_rank, 2L)
  expect_setequal(sub("^.*[.]", "", fit$curve_info$condition), c("A", "B"))
  for (curve_index in seq_along(fit$curves)) {
    condition <- sub("^.*[.]", "", fit$curve_info$condition[[curve_index]])
    for (voxel in seq_along(fit$voxels)) {
      truth <- truth_shape * amplitudes[condition, voxel]
      expect_gt(stats::cor(fit$estimate[, curve_index, voxel], truth), 0.9)
    }
  }
})

test_that("estimate_hrf fails clearly on invalid numerical contracts", {
  fixture <- .make_estimate_hrf_fixture(n_voxels = 1L, include_fixed = FALSE)

  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~run, dataset = fixture$dataset,
      lambda = "gcv", lambda_grid = c(0, -1)
    ),
    "non-negative"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), fixed = 1, block = ~run,
      dataset = fixture$dataset
    ),
    "fixed must"
  )

  nonfinite_data <- fixture$Y
  nonfinite_data[3, 1] <- NA_real_
  nonfinite_dataset <- fmridataset::matrix_dataset(
    nonfinite_data,
    TR = 1,
    run_length = nrow(nonfinite_data),
    event_table = fixture$events
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~run, dataset = nonfinite_dataset
    ),
    "non-finite"
  )

  aliased_events <- data.frame(
    onset = seq(10, 80, by = 10),
    condition = factor(rep(c("A", "B"), 4L)),
    condition_copy = factor(rep(c("A", "B"), 4L)),
    run = 1L
  )
  aliased_dataset <- fmridataset::matrix_dataset(
    matrix(rnorm(120), ncol = 1L),
    TR = 1,
    run_length = 120L,
    event_table = aliased_events
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition) + hrf(condition_copy),
      block = ~run, dataset = aliased_dataset,
      rsam = 0:20, k = 8L
    ),
    "rank deficient"
  )
})
