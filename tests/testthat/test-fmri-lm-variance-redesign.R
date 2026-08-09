test_that("sandwich covariance matches the HC1 oracle", {
  skip_if_not_installed("sandwich")
  set.seed(8101)
  n <- 90L
  X <- cbind(1, x = rnorm(n), z = runif(n, -1, 1))
  y <- drop(X %*% c(0.5, -0.25, 0.75) + rnorm(n) * (0.5 + abs(X[, 2])))
  fit <- stats::lm.fit(X, y)

  got <- fmrireg:::.fmri_lm_score_covariance(
    X, fit$residuals, segments = list(seq_len(n)),
    max_lag = 0L, taper = "none", debias = TRUE
  )
  lm_object <- stats::lm(y ~ x + z, data = data.frame(y, x = X[, 2], z = X[, 3]))
  expected <- sandwich::vcovHC(lm_object, type = "HC1")

  expect_equal(got, unname(expected), tolerance = 1e-8)
  expect_true(all(eigen(got, symmetric = TRUE, only.values = TRUE)$values >= -1e-10))
})

test_that("HAC meat preserves run boundaries", {
  set.seed(8102)
  n <- 24L
  X <- cbind(1, scale(seq_len(n)))
  residual <- rnorm(n)
  segments <- list(1:12, 13:24)
  H <- 2L

  got <- fmrireg:::.fmri_lm_score_covariance(
    X, residual, segments, max_lag = H,
    taper = "none", debias = FALSE
  )

  score <- X * residual
  meat <- crossprod(score)
  for (idx in segments) {
    S <- score[idx, , drop = FALSE]
    for (h in seq_len(H)) {
      gamma <- crossprod(S[(h + 1L):nrow(S), , drop = FALSE],
                         S[seq_len(nrow(S) - h), , drop = FALSE])
      meat <- meat + gamma + t(gamma)
    }
  }
  meat <- fmrireg:::.fmri_lm_psd(meat)
  bread <- solve(crossprod(X))
  expected <- fmrireg:::.fmri_lm_psd(bread %*% meat %*% bread)

  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("HAC Satterthwaite df reflects covariance-estimator uncertainty", {
  n <- 40L
  X <- matrix(1, nrow = n, ncol = 1L)
  residual <- rep(c(-1, 1), n / 2L)
  segments <- list(seq_len(n))

  iid_df <- fmrireg:::.fmri_lm_hac_df(
    X, residual, segments, max_lag = 0L, taper = "none"
  )
  lagged_df <- fmrireg:::.fmri_lm_hac_df(
    X, residual, segments, max_lag = 4L, taper = "none"
  )

  expect_equal(iid_df, n - 1L, tolerance = 1e-8)
  expect_true(is.finite(lagged_df))
  expect_gte(lagged_df, 1)
  expect_lte(lagged_df, n - 1L)
})

test_that("end-to-end HAC regenerates statistics with voxelwise df", {
  set.seed(8103)
  n <- 100L
  events <- data.frame(
    onset = seq(5, 75, by = 10),
    condition = factor(rep(c("A", "B"), 4)),
    run = 1L
  )
  Y <- replicate(3, as.numeric(arima.sim(list(ar = 0.45), n = n)))
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    variance = variance_spec("hac", max_lag = 4L, taper = "tukey",
                             df = "satterthwaite")
  )
  fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 2L)
  )

  vm <- variance_model(fit)
  beta <- fit$result$betas$data[[1L]]
  stat <- beta$stat[[1L]]
  prob <- beta$prob[[1L]]

  expect_identical(vm$method, "hac")
  expect_identical(vm$covariance_scope, "voxel")
  expect_length(vm$covariance, ncol(Y))
  expect_true(all(vapply(vm$covariance, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(is.finite(fit$result$df$inference)))
  expect_true(any(abs(fit$result$df$inference - fit$result$df$nominal) > 1e-6))
  for (v in seq_len(ncol(Y))) {
    expected <- 2 * stats::pt(-abs(stat[v, ]), df = fit$result$df$inference[v])
    expect_equal(prob[v, ], expected, tolerance = 1e-12)
  }

  weights <- structure(c(1, -1), colind = fit$result$event_indices)
  posthoc <- fit_contrasts(fit, list(A_vs_B = weights))$A_vs_B
  event_idx <- fit$result$event_indices
  manual_se <- vapply(seq_len(ncol(Y)), function(v) {
    cv <- vm$covariance[[v]][event_idx, event_idx, drop = FALSE]
    sqrt(drop(weights %*% cv %*% weights))
  }, numeric(1))
  expect_equal(posthoc$se, manual_se, tolerance = 1e-10)
  expect_equal(posthoc$df.inference, fit$result$df$inference)

  tidy_beta <- tidy(fit, type = "estimates")
  expect_true("df_inference" %in% names(tidy_beta))
  expect_equal(unique(tidy_beta$df_inference), unique(fit$result$df$inference))
  expect_equal(
    fmrireg:::.extract_degrees_of_freedom(fit),
    fit$result$df$inference
  )
  metadata <- fmrireg:::.fmri_lm_bids_inference_metadata(fit)
  expect_identical(metadata$VarianceMethod, "hac")
  expect_equal(metadata$InferenceDegreesOfFreedom$Min,
               min(fit$result$df$inference))
})

test_that("robust sandwich retains weights and uses their df adjustment", {
  set.seed(8104)
  n <- 80L
  events <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  Y <- matrix(rnorm(n * 2), n, 2)
  Y[c(12, 39, 61), ] <- Y[c(12, 39, 61), ] + 15
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    robust = robust_spec("huber", max_iter = 5L),
    variance = variance_spec("sandwich")
  )
  fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 1L)
  )

  expect_length(fit$result$fit_state$robust_weights, n)
  expect_true(any(fit$result$fit_state$robust_weights < 0.5))
  expect_true(all(fit$result$df$inference < fit$result$df$nominal))
  expect_true(all(is.finite(unlist(standard_error(fit)))))
})

test_that("AR fitting feeds post-whitening residuals into Satterthwaite inference", {
  set.seed(8105)
  n <- 90L
  events <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55, 65, 75),
    condition = factor(rep(c("A", "B"), 4)),
    run = 1L
  )
  Y <- replicate(2, as.numeric(arima.sim(list(ar = 0.6), n = n)))
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    noise = noise_spec("ar1", iter_gls = 2L),
    variance = variance_spec("model", max_lag = 4L,
                             df = "satterthwaite")
  )
  fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 1L)
  )

  expect_length(fit$result$fit_state$ar_parameters, 1L)
  expect_true(all(is.finite(unlist(fit$result$fit_state$ar_parameters))))
  expect_true(all(is.finite(fit$result$df$inference)))
  expect_identical(variance_model(fit)$metadata$noise$struct, "ar1")
  expect_equal(
    fit$result$betas$data[[1L]]$prob[[1L]],
    2 * stats::pt(
      -abs(fit$result$betas$data[[1L]]$stat[[1L]]),
      df = fit$result$df$inference
    ),
    tolerance = 1e-12
  )
})

test_that("runwise adaptive covariance is rejected until run-level combination exists", {
  dset <- .demo_matrix_dataset()
  control <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    variance = variance_spec("hac", max_lag = 2L)
  )
  expect_error(
    fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset,
            control = control, compute = compute_spec()),
    "require.*scope = 'joint'"
  )
})
