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
  dset <- matrix_frame(matrix(0, 40, 1), TR = 1, run_length = c(20L, 20L))
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
