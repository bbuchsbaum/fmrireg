# strategies preprocess NA/soft-subspace edges + meta tau2 singular/design edges.

test_that("preprocess_run_data errors on non-finite Y and missing soft nuisance", {
  set.seed(31)
  n <- 40L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * V), n, V)

  cfg <- fmri_lm_control()
  Y_bad <- Y
  Y_bad[1, 1] <- NA_real_
  expect_error(
    fmrireg:::preprocess_run_data(X, Y_bad, cfg, run_num = 1L, dataset = NULL),
    "Non-finite values"
  )

  cfg_prop <- fmri_lm_control(na_action = "propagate")
  out_prop <- fmrireg:::preprocess_run_data(X, Y_bad, cfg_prop, run_num = 1L)
  expect_equal(dim(out_prop$Y), dim(Y_bad))

  cfg_soft <- fmri_lm_control()
  cfg_soft$soft_subspace$enabled <- TRUE
  expect_error(
    fmrireg:::preprocess_run_data(X, Y, cfg_soft),
    "nuisance_matrix or nuisance_mask"
  )

  # Soft subspace with provided nuisance + run subsetting
  sframe <- fmrihrf::sampling_frame(c(20L, 20L), TR = 1)
  dset <- list(sampling_frame = sframe)
  N <- cbind(rnorm(40), rnorm(40))
  cfg_soft$soft_subspace$nuisance_matrix <- N
  cfg_soft$soft_subspace$lambda <- 0.25
  out_run <- fmrireg:::preprocess_run_data(
    X[1:20, , drop = FALSE], Y[1:20, , drop = FALSE],
    cfg_soft, run_num = 1L, dataset = dset
  )
  expect_equal(nrow(out_run$X), 20L)

  # Volume weights as matrix with run columns
  cfg_w <- fmri_lm_control()
  cfg_w$volume_weights$enabled <- TRUE
  cfg_w$volume_weights$weights <- cbind(runif(20, 0.5, 1), runif(20, 0.5, 1))
  out_w <- fmrireg:::preprocess_run_data(
    X[1:20, , drop = FALSE], Y[1:20, , drop = FALSE],
    cfg_w, run_num = 2L, dataset = dset
  )
  expect_equal(length(out_w$preprocess_info$volume_weights), 20L)
})

test_that("estimate_tau2 returns 0 for singular designs and covers reml", {
  # Singular design (duplicate columns)
  y <- rnorm(6)
  se <- rep(0.2, 6)
  X_sing <- cbind(rep(1, 6), rep(1, 6))
  expect_equal(fmrireg:::estimate_tau2(y, se, X_sing, method = "dl"), 0)

  set.seed(32)
  y2 <- rnorm(12, 0.3, 0.5)
  se2 <- runif(12, 0.1, 0.3)
  X2 <- cbind(1, rep(0:1, each = 6))
  tau_reml <- fmrireg:::estimate_tau2(y2, se2, X2, method = "reml")
  expect_true(is.finite(tau_reml) && tau_reml >= 0)

  # Homogeneous -> dl zero branch
  y0 <- rep(0.1, 8)
  se0 <- rep(0.2, 8)
  X0 <- cbind(1, rep(0:1, each = 4))
  expect_equal(fmrireg:::estimate_tau2(y0, se0, X0, method = "dl"), 0)
})
