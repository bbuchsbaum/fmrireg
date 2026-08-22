test_that("built-in joint fitting rejects voxelwise temporal covariance", {
  n <- 60L
  events <- data.frame(
    onset = c(10, 30), condition = factor(c("A", "B")), run = 1L
  )
  dset <- fmridataset::matrix_dataset(
    matrix(rnorm(2L * n), nrow = n),
    TR = 1, run_length = n, event_table = events
  )
  control <- fmri_lm_control(noise = noise_spec("ar1", voxelwise = TRUE))

  expect_true(control$noise$voxelwise)
  expect_identical(control$estimation$scope, "joint")
  expect_error(
    fmri_lm(
      onset ~ hrf(condition), block = ~run, dataset = dset,
      control = control, compute = compute_spec(backend = "matrix")
    ),
    "Built-in voxelwise temporal covariance currently requires.*runwise_meta"
  )

  runwise_control <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("ar1", voxelwise = TRUE)
  )
  expect_invisible(fmrireg:::.validate_builtin_voxelwise_noise(runwise_control))
})

test_that("built-in voxelwise noise rejects unsupported adjacent options", {
  voxelwise_control <- function(noise, weights = NULL, projection = NULL) {
    fmri_lm_control(
      estimation = estimation_spec("runwise_meta"),
      noise = noise,
      weights = weights,
      projection = projection
    )
  }

  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(noise_spec("iid", voxelwise = TRUE))
    ),
    "requires an AR structure"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(noise_spec("ar1", pooling = "global", voxelwise = TRUE))
    ),
    "requires `pooling = 'run'`"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(noise_spec("ar1", voxelwise = TRUE, censor = 5L))
    ),
    "does not yet support censoring"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(noise_spec("ar1", voxelwise = TRUE, iter_gls = 2L))
    ),
    "requires `iter_gls = 1`"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      fmri_lm_control(
        estimation = estimation_spec("runwise_meta"),
        noise = noise_spec("ar1", voxelwise = TRUE),
        robust = robust_spec("huber", reestimate_phi = TRUE)
      )
    ),
    "does not yet support robust AR re-estimation"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(
        noise_spec("ar1", voxelwise = TRUE),
        weights = weights_spec("tukey")
      )
    ),
    "does not yet support volume weighting"
  )
  expect_error(
    fmrireg:::.validate_builtin_voxelwise_noise(
      voxelwise_control(
        noise_spec("ar1", voxelwise = TRUE),
        projection = projection_spec(
          "soft_subspace", nuisance_matrix = matrix(rnorm(120), ncol = 2L)
        )
      )
    ),
    "soft-subspace projection"
  )
})

test_that("shared and voxelwise runwise estimators match their stated oracles", {
  set.seed(20260822)
  n <- 180L
  phi <- c(0.05, 0.25, 0.45, 0.65)
  Y <- vapply(
    phi,
    function(value) as.numeric(stats::arima.sim(list(ar = value), n = n)),
    numeric(n)
  )
  events <- data.frame(
    onset = seq(10, 150, by = 20),
    condition = factor(rep(c("A", "B"), 4L)),
    run = 1L
  )
  dset <- fmridataset::matrix_dataset(
    Y, TR = 1, run_length = n, event_table = events
  )
  control_for <- function(voxelwise, shared_estimator = "pooled_acvf") {
    fmri_lm_control(
      estimation = estimation_spec("runwise_meta"),
      noise = noise_spec(
        "ar1", voxelwise = voxelwise,
        shared_estimator = shared_estimator
      ),
      variance = variance_spec("model")
    )
  }

  shared_fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control_for(FALSE), compute = compute_spec(backend = "matrix")
  )
  mean_series_fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control_for(FALSE, "mean_series"),
    compute = compute_spec(backend = "matrix")
  )
  voxelwise_fit <- fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = dset,
    control = control_for(TRUE), compute = compute_spec(backend = "matrix")
  )

  tmats <- fmridesign::term_matrices(voxelwise_fit$model, 1L)
  data_env <- list2env(tmats)
  data_env$.y <- Y
  X <- stats::model.matrix(fmrireg:::get_formula(voxelwise_fit$model), data_env)
  proj <- fmrireg:::.fast_preproject(X)
  ols <- fmrireg:::solve_glm_core(
    fmrireg:::glm_context(X = X, Y = Y, proj = proj),
    return_fitted = TRUE
  )
  residuals <- Y - ols$fitted
  shared_oracle <- fmrireg:::.estimate_ar_parameters_routed(
    residuals, 1L, design = X
  )
  mean_series_oracle <- fmrireg:::.estimate_ar_parameters_routed(
    rowMeans(residuals), 1L, design = X
  )
  voxelwise_oracle <- vapply(
    seq_len(ncol(Y)),
    function(v) fmrireg:::estimate_ar_parameters(residuals[, v], 1L),
    numeric(1)
  )

  shared <- unlist(ar_parameters(shared_fit, scope = "raw"), use.names = FALSE)
  mean_series <- unlist(
    ar_parameters(mean_series_fit, scope = "raw"), use.names = FALSE
  )
  voxelwise <- unlist(ar_parameters(voxelwise_fit, scope = "raw"), use.names = FALSE)
  expect_length(shared, 1L)
  expect_length(mean_series, 1L)
  expect_length(voxelwise, ncol(Y))
  expect_equal(shared, shared_oracle, tolerance = 1e-12)
  expect_equal(mean_series, mean_series_oracle, tolerance = 1e-12)
  expect_equal(voxelwise, voxelwise_oracle, tolerance = 1e-12)
  expect_true(attr(voxelwise_fit, "executed_config")$noise$voxelwise)
  expect_identical(
    attr(shared_fit, "executed_config")$noise$shared_estimator,
    "pooled_acvf"
  )
})
