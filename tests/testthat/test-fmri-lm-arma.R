test_that("MA specifications expose a narrow supported contract", {
  arma <- noise_spec("ar1", q = 1L, iter_gls = 2L)
  expect_identical(arma$p, 1L)
  expect_identical(arma$q, 1L)

  ma <- noise_spec("iid", q = 2L)
  expect_null(ma$p)
  expect_identical(ma$q, 2L)

  expect_error(noise_spec("ar1", q = 1L, pooling = "global"),
               "MA terms currently require.*pooling = 'run'")
  expect_error(noise_spec("ar1", q = 1L, pooling = "parcel"),
               "MA terms currently require.*pooling = 'run'")
  expect_error(noise_spec("ar1", q = 1L, voxelwise = TRUE),
               "MA terms do not yet support voxelwise")
  expect_error(noise_spec("ar1", q = 1L, censor = 10L),
               "MA terms do not yet support censoring")
  expect_error(noise_spec("ar1", q = 1L, iter_gls = 0L),
               "MA terms require `iter_gls >= 1`")
  expect_error(
    fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec("ar1", q = 1L)
    ),
    "MA terms currently require.*runwise_meta"
  )
  expect_error(
    fmri_lm_control(
      estimation = estimation_spec("runwise_meta"),
      noise = noise_spec("ar1", q = 1L),
      robust = robust_spec("huber")
    ),
    "MA terms cannot yet be combined with robust"
  )
})

test_that("fmrireg maps q > 0 to fmriAR ARMA estimation", {
  cfg <- fmrireg:::.fmrireg_to_fmriAR_config(
    list(struct = "ar1", p = 1L, q = 2L, iter_gls = 3L),
    n_runs = 1L
  )

  expect_identical(cfg$method, "arma")
  expect_identical(cfg$p, 1L)
  expect_identical(cfg$q, 2L)
  expect_identical(cfg$iter, 3L)
})

test_that("iterative ARMA GLS carries theta through an fmriAR whitening oracle", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(82201)
  n <- 500L
  n_vox <- 3L
  X <- cbind(intercept = 1, trend = as.numeric(scale(seq_len(n))))
  beta <- matrix(c(0.5, 0.25), nrow = 2L, ncol = n_vox)
  errors <- replicate(
    n_vox,
    as.numeric(stats::arima.sim(list(ar = 0.45, ma = 0.35), n = n))
  )
  Y <- X %*% beta + errors

  out <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X = X,
    Y = Y,
    cfg = list(
      struct = "ar1", p = 1L, q = 1L, iter_gls = 2L,
      pooling = "run", exact_first = FALSE
    ),
    max_iter = 2L
  )
  oracle <- fmriAR::whiten_apply(out$plan, X = X, Y = Y, parallel = FALSE)

  expect_s3_class(out$plan, "fmriAR_plan")
  expect_identical(unname(out$plan$order[c("p", "q")]), c(1L, 1L))
  expect_length(unlist(out$ar_coef), 1L)
  expect_length(unlist(out$ma_coef), 1L)
  expect_true(all(is.finite(unlist(out$ar_coef))))
  expect_true(all(is.finite(unlist(out$ma_coef))))
  expect_equal(out$X_white, oracle$X, tolerance = 1e-12)
  expect_equal(out$Y_white, oracle$Y, tolerance = 1e-12)

  beta_gls <- fmrireg:::.fast_preproject(out$X_white)$Pinv %*% out$Y_white
  whitened_resid <- out$Y_white - out$X_white %*% beta_gls
  raw_resid <- Y - X %*% base::qr.solve(X, Y)
  lag1 <- function(z) stats::acf(rowMeans(z), plot = FALSE, lag.max = 1)$acf[2L]
  expect_lt(abs(lag1(whitened_resid)), abs(lag1(raw_resid)))
})

test_that("pure MA fitting is not mistaken for IID", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(82202)
  n <- 350L
  X <- cbind(1, as.numeric(scale(seq_len(n))))
  Y <- matrix(
    X %*% c(0.5, -0.2) +
      as.numeric(stats::arima.sim(list(ma = 0.5), n = n)),
    ncol = 1L
  )

  out <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X = X,
    Y = Y,
    cfg = list(
      struct = "iid", p = NULL, q = 1L, iter_gls = 1L,
      pooling = "run", exact_first = FALSE
    )
  )

  expect_s3_class(out$plan, "fmriAR_plan")
  expect_identical(unname(out$plan$order[["p"]]), 0L)
  expect_identical(unname(out$plan$order[["q"]]), 1L)
  expect_length(unlist(out$ma_coef), 1L)
  expect_false(isTRUE(all.equal(out$X_white, X)))
})

test_that("ARMA effective degrees of freedom include theta", {
  n <- 240L
  n_beta <- 4L
  phi <- 0.4
  theta <- 0.3
  rho <- stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1L)[-1L]
  lag <- seq_along(rho)
  expected <- n / (1 + 2 * sum((1 - lag / n) * rho)) - n_beta

  actual <- fmrireg:::.compute_ar_effective_df_compat(
    n, n_beta, ar_coef = list(phi), ma_coef = list(theta)
  )
  ar_only <- fmrireg:::.compute_ar_effective_df_compat(
    n, n_beta, ar_coef = list(phi)
  )

  expect_equal(actual, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(actual, ar_only)))
})

test_that("public runwise ARMA fitting returns per-run parameters and inference", {
  skip_if_not_installed("fmriAR", minimum_version = "0.3.3")

  set.seed(82203)
  n <- 180L
  n_run <- 2L
  events <- do.call(rbind, lapply(seq_len(n_run), function(run) {
    data.frame(
      onset = seq(8, 152, by = 16),
      condition = factor(rep(c("A", "B"), 5L), levels = c("A", "B")),
      run = run
    )
  }))
  Y <- do.call(rbind, lapply(seq_len(n_run), function(run) {
    replicate(
      2L,
      as.numeric(stats::arima.sim(list(ar = 0.4, ma = 0.3), n = n))
    )
  }))
  dset <- fmridataset::matrix_dataset(
    Y, TR = 1, run_length = rep(n, n_run), event_table = events
  )
  control <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("ar1", q = 1L, iter_gls = 2L),
    variance = variance_spec("model")
  )

  fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control, compute = compute_spec(voxel_chunks = 1L)
  )

  expect_s3_class(fit, "fmri_lm")
  expect_length(unlist(ar_parameters(fit, scope = "raw")), n_run)
  expect_length(unlist(ma_parameters(fit, scope = "raw")), n_run)
  expect_length(ma_parameters(fit, scope = "per_run"), n_run)
  expect_length(ma_parameters(fit, scope = "average"), 1L)
  expect_true(all(is.finite(unlist(ma_parameters(fit, scope = "raw")))))
  expect_length(fit$result$fit_state$temporal_diagnostics, n_run)
  expect_true(all(vapply(
    fit$result$fit_state$temporal_diagnostics,
    function(x) x$iterations >= 1L && x$iterations <= 2L,
    logical(1)
  )))
  expect_true(all(is.finite(unlist(standard_error(fit)))))
  expect_identical(fit$result$fit_state$ma_parameters, fit$result$ma_coef)
  expect_identical(variance_model(fit)$metadata$noise$q, 1L)
})

test_that("unsupported fitters reject MA terms before dispatch", {
  control <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("ar1", q = 1L)
  )

  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "ar_only", control, capabilities = list(ma = FALSE)
    ),
    "ar_only does not support MA terms"
  )
})
