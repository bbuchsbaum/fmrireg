# Independent contract tests for shared AR estimation.  These deliberately use
# small algebraic oracles rather than reproducing the fmriAR implementation.

.manual_pooled_ar1 <- function(residuals) {
  residuals <- as.matrix(residuals)
  residuals <- sweep(residuals, 2L, colMeans(residuals), "-")
  n <- nrow(residuals)
  gamma0 <- sum(residuals * residuals) / (ncol(residuals) * n)
  gamma1 <- sum(
    residuals[-1L, , drop = FALSE] *
      residuals[-n, , drop = FALSE]
  ) / (ncol(residuals) * (n - 1L))
  max(-0.99, min(0.99, gamma1 / gamma0))
}

.oracle_bias_map <- function(design, max_lag) {
  design <- as.matrix(design)
  n <- nrow(design)
  qx <- qr(design)
  Q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
  M <- diag(n) - tcrossprod(Q)
  A <- matrix(0, max_lag + 1L, max_lag + 1L)

  for (k in 0:max_lag) {
    S <- matrix(0, n, n)
    if (k == 0L) {
      diag(S) <- 1
    } else {
      lo <- seq_len(n - k)
      S[cbind(lo, lo + k)] <- 1
      S[cbind(lo + k, lo)] <- 1
    }
    expected_residual_product <- M %*% S %*% M
    for (h in 0:max_lag) {
      lo <- seq_len(n - h)
      A[h + 1L, k + 1L] <- mean(
        expected_residual_product[cbind(lo, lo + h)]
      )
    }
  }
  A
}

.oracle_expected_corrected_ar1 <- function(design, phi, max_lag) {
  n <- nrow(design)
  qx <- qr(design)
  Q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
  M <- diag(n) - tcrossprod(Q)
  Sigma <- phi^abs(outer(seq_len(n), seq_len(n), "-"))
  residual_covariance <- M %*% Sigma %*% M
  raw <- vapply(0:max_lag, function(h) {
    lo <- seq_len(n - h)
    mean(residual_covariance[cbind(lo, lo + h)])
  }, numeric(1L))
  corrected <- solve(.oracle_bias_map(design, max_lag), raw)
  corrected[2L] / corrected[1L]
}

test_that("fixed-order pooled and mean-series estimands remain distinct", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(790)
  n <- 400L
  common <- as.numeric(arima.sim(list(ar = 0.75), n = n))
  opposing <- 3 * as.numeric(arima.sim(list(ar = 0.05), n = n))
  residuals <- cbind(common + opposing, common - opposing)

  pooled_cfg <- list(struct = "ar1", shared_estimator = "pooled_acvf")
  mean_cfg <- list(struct = "ar1", shared_estimator = "mean_series")
  pooled <- .estimate_shared_ar_parameters(residuals, 1L, pooled_cfg)
  mean_series <- .estimate_shared_ar_parameters(residuals, 1L, mean_cfg)

  expect_equal(pooled, .manual_pooled_ar1(residuals), tolerance = 1e-12)
  expect_equal(
    mean_series,
    .manual_pooled_ar1(matrix(rowMeans(residuals), ncol = 1L)),
    tolerance = 1e-12
  )
  expect_gt(abs(mean_series - pooled), 0.2)
})

test_that("shared run pooling preserves every run-specific AR vector", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(795)
  n_run <- 140L
  run_indices <- list(seq_len(n_run), n_run + seq_len(n_run))
  residuals <- rbind(
    replicate(4L, as.numeric(arima.sim(list(ar = 0.15), n = n_run))),
    replicate(4L, as.numeric(arima.sim(list(ar = 0.7), n = n_run)))
  )
  cfg <- list(struct = "ar1", global = FALSE,
              shared_estimator = "pooled_acvf")

  actual <- .estimate_shared_ar_parameters(
    residuals, 1L, cfg, run_indices = run_indices
  )
  expected <- fmriAR::fit_noise(
    residuals, runs = rep(1:2, each = n_run),
    method = "ar", p = 1L, p_max = 1L, pooling = "run"
  )$phi

  expect_type(actual, "list")
  expect_length(actual, 2L)
  expect_equal(actual, lapply(expected, as.numeric), tolerance = 1e-12)
  expect_gt(abs(actual[[2L]] - actual[[1L]]), 0.3)
})

test_that("mean-series pooling retains the matching OLS correction", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(2024)
  n <- 235L
  X <- cbind(1, poly(seq_len(n), 4), matrix(rnorm(n * 20L), n))
  Y <- matrix(rnorm(n * 12L), n)
  residuals <- qr.resid(qr(X), Y)
  mean_cfg <- list(struct = "ar1", shared_estimator = "mean_series")

  mean_after_projection <- rowMeans(residuals)
  projection_after_mean <- drop(qr.resid(qr(X), rowMeans(Y)))
  expect_equal(mean_after_projection, projection_after_mean, tolerance = 1e-12)
  expect_identical(.shared_ar_correction_design(residuals, mean_cfg, X), X)

  budget <- .ar_correction_lag_budget(X, target_order = 1L)
  expected <- fmriAR::fit_noise(
    matrix(mean_after_projection, ncol = 1L),
    method = "ar", p = 1L, p_max = 1L, pooling = "global",
    design = X, correction_max_lag = budget
  )
  actual <- .estimate_shared_ar_parameters(
    residuals, 1L, mean_cfg, design = X
  )
  expect_equal(actual, as.numeric(expected$phi[[1L]]), tolerance = 1e-12)
})

test_that("AR correction budget covers projection-induced covariance tails", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(124)
  n <- 120L
  X <- cbind(1, poly(seq_len(n), 4), matrix(rnorm(n * 34L), n))
  budget <- .ar_correction_lag_budget(X, target_order = 1L)

  expect_identical(budget, 5L)
  phi <- 0.8
  at_order <- .oracle_expected_corrected_ar1(X, phi, max_lag = 1L)
  with_tail <- .oracle_expected_corrected_ar1(X, phi, max_lag = budget)
  expect_lt(abs(with_tail - phi), abs(at_order - phi))
  expect_lt(abs(with_tail - phi), 0.03)

  # The policy is bounded for cost, never requests less than the fitted order,
  # and contracts safely when a short, rich design has little residual df.
  X_short <- cbind(1, poly(seq_len(18L), 5), diag(18L)[, 1:7, drop = FALSE])
  short_budget <- .ar_correction_lag_budget(X_short, target_order = 2L)
  expect_gte(short_budget, 2L)
  expect_lte(short_budget, 25L)
})

test_that("block-diagonal design exactly represents separate runwise OLS fits", {
  set.seed(791)
  n1 <- 45L
  n2 <- 37L
  X1 <- cbind(1, scale(seq_len(n1)), rnorm(n1))
  X2 <- cbind(1, scale(seq_len(n2)), rnorm(n2))
  Y1 <- matrix(rnorm(n1 * 4L), n1)
  Y2 <- matrix(rnorm(n2 * 4L), n2)
  separate_residuals <- rbind(qr.resid(qr(X1), Y1), qr.resid(qr(X2), Y2))

  block_design <- .block_diagonal_design(list(X1, X2))
  block_residuals <- qr.resid(qr(block_design), rbind(Y1, Y2))
  shared_design_residuals <- qr.resid(qr(rbind(X1, X2)), rbind(Y1, Y2))

  expect_equal(dim(block_design), c(n1 + n2, ncol(X1) + ncol(X2)))
  expect_true(all(block_design[
    seq_len(n1), ncol(X1) + seq_len(ncol(X2)), drop = FALSE
  ] == 0))
  expect_true(all(block_design[
    n1 + seq_len(n2), seq_len(ncol(X1)), drop = FALSE
  ] == 0))
  expect_equal(block_residuals, separate_residuals, tolerance = 1e-12)
  expect_gt(max(abs(shared_design_residuals - separate_residuals)), 0.05)

  # Different per-run model widths are supported without recycling or padding
  # columns into the wrong run's coefficient space.
  X2_wide <- cbind(X2, rnorm(n2), X2[, 2L])
  heterogeneous <- .block_diagonal_design(list(X1, X2_wide))
  expect_equal(dim(heterogeneous), c(n1 + n2, ncol(X1) + ncol(X2_wide)))
})

test_that("production global AR routes the block residual operator", {
  dset <- .demo_matrix_dataset()
  model <- .demo_fmri_model()
  cfg <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("ar1", pooling = "global")
  )
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .estimate_shared_ar_parameters = function(residuals, ar_order, ar_opts,
                                              run_indices = NULL, censor = NULL,
                                              design = NULL) {
      captured$residuals <- residuals
      captured$design <- design
      captured$run_indices <- run_indices
      stop("captured global AR", call. = FALSE)
    },
    .package = "fmrireg"
  )

  expect_error(
    suppressWarnings(fmri_lm_fit(
      model, dset, strategy = "runwise", cfg = cfg,
      use_fast_path = TRUE, progress = FALSE
    )),
    "captured global AR"
  )

  expect_identical(unname(lengths(captured$run_indices)), c(4L, 4L))
  expect_equal(dim(captured$design), c(8L, 12L))
  first <- captured$run_indices[[1L]]
  second <- captured$run_indices[[2L]]
  expect_true(all(captured$design[first, 7:12, drop = FALSE] == 0))
  expect_true(all(captured$design[second, 1:6, drop = FALSE] == 0))
  expect_lt(max(abs(crossprod(captured$design, captured$residuals))), 1e-10)
})

test_that("AR whitening resets exactly at every run boundary", {
  skip_if_not_installed("fmriAR")

  set.seed(792)
  lengths <- c(13L, 11L)
  run_indices <- list(seq_len(lengths[1L]), lengths[1L] + seq_len(lengths[2L]))
  X <- cbind(1, rnorm(sum(lengths)))
  Y <- matrix(rnorm(sum(lengths) * 3L), ncol = 3L)

  combined <- ar_whiten_transform(
    X, Y, phi = 0.55, exact_first = TRUE, run_indices = run_indices
  )
  separate <- lapply(run_indices, function(idx) {
    ar_whiten_transform(X[idx, , drop = FALSE], Y[idx, , drop = FALSE],
                        phi = 0.55, exact_first = TRUE)
  })
  expect_equal(combined$X, do.call(rbind, lapply(separate, `[[`, "X")),
               tolerance = 1e-12)
  expect_equal(combined$Y, do.call(rbind, lapply(separate, `[[`, "Y")),
               tolerance = 1e-12)

  crossing <- ar_whiten_transform(X, Y, phi = 0.55, exact_first = TRUE)
  first_second_run <- lengths[1L] + 1L
  expect_gt(max(abs(crossing$Y[first_second_run, ] -
                    combined$Y[first_second_run, ])), 0.05)
})

test_that("run-index labels fail closed unless they form an exact partition", {
  expect_identical(
    .build_run_labels(6L, list(1:3, 4:6)),
    c(1L, 1L, 1L, 2L, 2L, 2L)
  )
  expect_error(
    .build_run_labels(6L, list(c(1L, 3L, 5L), c(2L, 4L, 6L))),
    "contiguous runs"
  )
  expect_error(.build_run_labels(6L, list(1:3, 3:6)), "exactly once")
  expect_error(.build_run_labels(6L, list(1:2, 4:6)), "every timepoint")
  expect_error(.build_run_labels(6L, list(0:2, 3:6)), "every timepoint")
  expect_error(.build_run_labels(6L, list(1:3, c(4, NA, 6))),
               "every timepoint")
})

test_that("low-rank global AR uses initial OLS once with full provenance", {
  set.seed(793)
  X <- cbind(1, rnorm(40L), rnorm(40L))
  Z <- matrix(rnorm(40L * 6L), 40L)
  run_indices <- list(1:17, 18:40)
  captured <- new.env(parent = emptyenv())
  captured$calls <- 0L

  local_mocked_bindings(
    .estimate_shared_ar_parameters = function(residuals, ar_order, ar_opts,
                                              run_indices = NULL, censor = NULL,
                                              design = NULL) {
      captured$calls <- captured$calls + 1L
      captured$residuals <- residuals
      captured$design <- design
      captured$run_indices <- run_indices
      0.4
    },
    ar_whiten_transform = function(X, Y, phi, exact_first = FALSE,
                                   censor = NULL, run_indices = NULL, ...) {
      captured$whitening_runs <- run_indices
      list(X = X, Y = Y)
    },
    .package = "fmrireg"
  )

  result <- .lowrank_whiten_initial_ols(
    X, Z, ar_order = 1L,
    ar_opts = list(struct = "ar1", iter_gls = 9L,
                   shared_estimator = "pooled_acvf"),
    run_indices = run_indices
  )

  expect_identical(captured$calls, 1L)
  expect_equal(captured$residuals, qr.resid(qr(X), Z), tolerance = 1e-12)
  expect_identical(captured$design, X)
  expect_identical(captured$run_indices, run_indices)
  expect_identical(captured$whitening_runs, run_indices)
  expect_equal(result$phi, 0.4)
})

test_that("low-rank parcel AR pools groups as columns, not fake time", {
  set.seed(794)
  X <- cbind(1, rnorm(30L))
  residuals <- list(a = rnorm(30L), b = rnorm(30L))
  sizes <- c(a = 10L, b = 30L)
  run_indices <- list(1:15, 16:30)
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()
  answers <- list(0.4, 0.2, 0.6)

  local_mocked_bindings(
    .estimate_shared_ar_parameters = function(residuals, ar_order, ar_opts,
                                              run_indices = NULL, censor = NULL,
                                              design = NULL) {
      i <- length(captured$calls) + 1L
      captured$calls[[i]] <- list(
        residuals = as.matrix(residuals), design = design,
        run_indices = run_indices
      )
      answers[[i]]
    },
    .package = "fmrireg"
  )

  result <- .lowrank_group_ar_estimates(
    residuals, sizes, ar_order = 1L,
    ar_opts = list(struct = "ar1", shared_estimator = "pooled_acvf"),
    shrink_c0 = 100L, design = X, run_indices = run_indices
  )

  expect_length(captured$calls, 3L)
  expect_equal(captured$calls[[1L]]$residuals,
               do.call(cbind, residuals), tolerance = 0)
  expect_true(all(vapply(captured$calls, function(x) identical(x$design, X),
                         logical(1L))))
  expect_true(all(vapply(captured$calls,
                         function(x) identical(x$run_indices, run_indices),
                         logical(1L))))
  expect_equal(result$global_phi, 0.4)
  expect_equal(result$phi_groups$a, (10 / 110) * 0.2 + (100 / 110) * 0.4)
  expect_equal(result$phi_groups$b, (30 / 130) * 0.6 + (100 / 130) * 0.4)
})

test_that("model run extraction is shared by specialized engines", {
  model <- .demo_fmri_model()
  runs <- .model_run_indices(model, n = 8L)
  expect_length(runs, 2L)
  expect_identical(unname(lengths(runs)), c(4L, 4L))
  expect_identical(unlist(runs, use.names = FALSE), seq_len(8L))
})
