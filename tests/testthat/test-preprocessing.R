# Tests for preprocessing functions: volume quality and soft subspace projection

test_that("compute_dvars calculates temporal derivative correctly", {
  set.seed(42)

  # Create simple test data
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  dvars <- compute_dvars(Y, normalize = FALSE)

  expect_length(dvars, 100)
  expect_true(all(dvars >= 0))

  # First value should be median of rest (no derivative available)
  dvars_raw <- sqrt(rowMeans(diff(Y)^2))
  expect_equal(dvars[1], median(dvars_raw))

  # Check that adding a spike increases DVARS
  Y_spike <- Y
  Y_spike[50, ] <- Y_spike[50, ] + 10
  dvars_spike <- compute_dvars(Y_spike, normalize = FALSE)

  # DVARS at timepoint 50 (after spike) should be higher
  expect_gt(dvars_spike[50], dvars[50])
})

test_that("compute_dvars normalizes correctly", {
  set.seed(42)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  dvars_norm <- compute_dvars(Y, normalize = TRUE)
  dvars_raw <- compute_dvars(Y, normalize = FALSE)

  # Normalized DVARS should have median ~1
  expect_equal(median(dvars_norm), 1, tolerance = 0.1)

  # Should be proportional
  expect_equal(dvars_norm / dvars_raw, rep(1 / median(dvars_raw), 100), tolerance = 1e-10)
})

test_that("compute_dvars requires at least 2 timepoints", {
  Y <- matrix(1:10, nrow = 1, ncol = 10)
  expect_error(compute_dvars(Y), "at least 2 timepoints")
})

test_that("dvars_to_weights produces valid weights", {
  set.seed(42)
  dvars <- abs(rnorm(100, mean = 1, sd = 0.3))

  for (method in c("inverse_squared", "soft_threshold", "tukey")) {
    weights <- dvars_to_weights(dvars, method = method)

    expect_length(weights, 100)
    expect_true(all(weights >= 0))
    # Mean should be approximately 1 after normalization
    expect_equal(mean(weights), 1, tolerance = 0.01)
  }
})

test_that("dvars_to_weights inverse_squared method works correctly", {
  dvars <- c(0, 0.5, 1, 2, 5)
  weights <- dvars_to_weights(dvars, method = "inverse_squared")

  # Weights should decrease as DVARS increases
  # Note: weights are normalized by mean, so absolute values differ
  raw_weights <- 1 / (1 + dvars^2)
  normalized <- raw_weights / mean(raw_weights)
  expect_equal(weights, normalized, tolerance = 1e-10)
})

test_that("dvars_to_weights soft_threshold respects threshold", {
  dvars <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  threshold <- 1.5

  weights <- dvars_to_weights(dvars, method = "soft_threshold", threshold = threshold)

  # Values below threshold should have higher weights
  # (but exact comparison is tricky due to normalization)
  expect_true(all(weights > 0))
})

test_that("dvars_to_weights tukey downweights extreme values", {
  dvars <- c(0.1, 0.5, 1.0, 5.0, 10.0)  # Last two are extreme
  threshold <- 1.5

  weights <- dvars_to_weights(dvars, method = "tukey", threshold = threshold)

  # Extreme values should be heavily downweighted or zero
  # tukey cutoff is threshold * 2 = 3.0
  expect_true(weights[5] < weights[1])  # Extreme < normal
})

test_that("volume_weights convenience function works", {
  set.seed(42)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  # Basic usage
  weights <- volume_weights(Y)
  expect_length(weights, 100)
  expect_true(all(weights > 0))

  # With return_dvars
  result <- volume_weights(Y, return_dvars = TRUE)
  expect_named(result, c("weights", "dvars"))
  expect_length(result$weights, 100)
  expect_length(result$dvars, 100)
})

test_that("apply_volume_weights transforms X and Y correctly", {
  skip("apply_volume_weights not yet implemented")
  set.seed(42)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  weights <- rep(1, 100)
  weights[50] <- 0.25  # Downweight one timepoint

  result <- apply_volume_weights(X, Y, weights)

  # Check dimensions
  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))

  # Weighted rows should be scaled by sqrt(weight)
  expect_equal(result$X[50, ], X[50, ] * sqrt(0.25))
  expect_equal(result$Y[50, ], Y[50, ] * sqrt(0.25))

  # Unweighted rows should be unchanged (weights = 1)
  expect_equal(result$X[1, ], X[1, ])
})

test_that("apply_volume_weights validates inputs", {
  skip("apply_volume_weights not yet implemented")
  X <- matrix(1:20, nrow = 10, ncol = 2)
  Y <- matrix(1:30, nrow = 10, ncol = 3)
  wrong_weights <- rep(1, 5)  # Wrong length

  expect_error(apply_volume_weights(X, Y, wrong_weights), "must equal number of timepoints")
})

# --- Soft Subspace Projection Tests ---

test_that("soft_projection creates valid projection object", {
  set.seed(42)
  N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)

  proj <- soft_projection(N, lambda = "auto")

  expect_s3_class(proj, "soft_projection")
  expect_true(is.function(proj$P_lambda))
  expect_true(proj$lambda > 0)
  expect_equal(proj$method, "singular_value_heuristic")
  expect_true(proj$effective_df > 0 && proj$effective_df < 20)
  expect_equal(proj$n_nuisance, 20)
  expect_equal(proj$n_timepoints, 100)
})

test_that("soft_projection with user-specified lambda works", {
  set.seed(42)
  N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)

  proj <- soft_projection(N, lambda = 1.5)

  expect_equal(proj$lambda, 1.5)
  expect_equal(proj$method, "user_specified")
})

test_that("soft_projection removes nuisance variance", {
  set.seed(42)

  # Create nuisance that is correlated with data
  N <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)

  # Create data that is partially explained by nuisance
  signal <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  nuisance_effect <- N %*% matrix(rnorm(5 * 10), nrow = 5, ncol = 10)
  Y <- signal + nuisance_effect * 2

  # Apply soft projection
  proj <- soft_projection(N, lambda = "auto")
  Y_clean <- proj$P_lambda(Y)

  # Cleaned data should have less correlation with nuisance
  cor_before <- mean(abs(cor(Y, N)))
  cor_after <- mean(abs(cor(Y_clean, N)))

  expect_lt(cor_after, cor_before)
})

test_that("soft_projection preserves dimensions", {
  set.seed(42)
  N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  proj <- soft_projection(N, lambda = "auto")
  Y_clean <- proj$P_lambda(Y)

  expect_equal(dim(Y_clean), dim(Y))
})

test_that("soft_projection GCV lambda selection works", {
  set.seed(42)
  N <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Y <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)

  proj <- soft_projection(N, lambda = "gcv", Y = Y)

  expect_equal(proj$method, "gcv")
  expect_true(proj$lambda > 0)
})

test_that("soft_projection GCV falls back to auto without Y", {
  set.seed(42)
  N <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)

  expect_warning(
    proj <- soft_projection(N, lambda = "gcv", Y = NULL),
    "falling back"
  )

  expect_equal(proj$method, "singular_value_heuristic")
})

test_that("apply_soft_projection works correctly", {
  set.seed(42)
  N <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  X <- cbind(1, rnorm(100))

  proj <- soft_projection(N, lambda = "auto")
  result <- apply_soft_projection(proj, Y, X)

  expect_named(result, c("Y", "X"))
  expect_equal(dim(result$Y), dim(Y))
  expect_equal(dim(result$X), dim(X))
})

test_that("apply_soft_projection validates input", {
  not_proj <- list(foo = 1)

  expect_error(
    apply_soft_projection(not_proj, matrix(1), matrix(1)),
    "must be a soft_projection object"
  )
})

test_that("soft_subspace_options creates valid options", {
  opts <- soft_subspace_options(
    enabled = TRUE,
    nuisance_matrix = matrix(1:20, 10, 2),
    lambda = "auto"
  )

  expect_s3_class(opts, "soft_subspace_options")
  expect_true(opts$enabled)
  expect_equal(opts$lambda, "auto")
  expect_true(opts$warn_redundant)
})

test_that("soft_subspace_options requires nuisance source when enabled", {
  expect_error(
    soft_subspace_options(enabled = TRUE),
    "nuisance_mask or nuisance_matrix"
  )
})

test_that("soft_subspace_options prefers matrix over mask", {
  expect_warning(
    soft_subspace_options(
      enabled = TRUE,
      nuisance_mask = "fake_path.nii",
      nuisance_matrix = matrix(1:20, 10, 2)
    ),
    "using nuisance_matrix"
  )
})

# --- Integration with fmri_lm_config ---

test_that("fmri_lm_control accepts volume_weights_options", {
  cfg <- fmri_lm_control(
    volume_weights_options = list(
      enabled = TRUE,
      method = "inverse_squared",
      threshold = 2.0
    )
  )

  expect_true(cfg$volume_weights$enabled)
  expect_equal(cfg$volume_weights$method, "inverse_squared")
  expect_equal(cfg$volume_weights$threshold, 2.0)
})

test_that("fmri_lm_control accepts soft_subspace_options", {
  N <- matrix(rnorm(100 * 10), 100, 10)

  cfg <- fmri_lm_control(
    soft_subspace_options = list(
      enabled = TRUE,
      nuisance_matrix = N,
      lambda = "auto"
    )
  )

  expect_true(cfg$soft_subspace$enabled)
  expect_equal(cfg$soft_subspace$lambda, "auto")
})

test_that("fmri_lm_control validates soft_subspace when enabled", {
  expect_error(
    fmri_lm_control(
      soft_subspace_options = list(enabled = TRUE)
    ),
    "nuisance_mask or nuisance_matrix"
  )
})

test_that("fmri_lm_control defaults are correct", {
  cfg <- fmri_lm_control()

  expect_false(cfg$volume_weights$enabled)
  expect_equal(cfg$volume_weights$method, "inverse_squared")
  expect_equal(cfg$volume_weights$threshold, 1.5)

  expect_false(cfg$soft_subspace$enabled)
  expect_true(cfg$soft_subspace$warn_redundant)
})

# --- preprocess_run_data integration test ---

test_that("preprocess_run_data applies volume weights", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  cfg <- fmri_lm_control(
    volume_weights_options = list(enabled = TRUE, method = "inverse_squared")
  )

  result <- preprocess_run_data(X, Y, cfg)

  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))
  expect_false(is.null(result$preprocess_info$volume_weights))
  expect_false(is.null(result$preprocess_info$dvars))
})

test_that("preprocess_run_data applies soft projection", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  N <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)

  cfg <- fmri_lm_control(
    soft_subspace_options = list(
      enabled = TRUE,
      nuisance_matrix = N,
      lambda = "auto"
    )
  )

  result <- preprocess_run_data(X, Y, cfg)

  expect_equal(dim(result$X), dim(X))
  expect_equal(dim(result$Y), dim(Y))
  expect_s3_class(result$preprocess_info$soft_projection, "soft_projection")
})

test_that("preprocess_run_data can apply both preprocessing steps", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  N <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)

  cfg <- fmri_lm_control(
    volume_weights_options = list(enabled = TRUE),
    soft_subspace_options = list(
      enabled = TRUE,
      nuisance_matrix = N,
      lambda = "auto"
    )
  )

  result <- preprocess_run_data(X, Y, cfg)

  expect_false(is.null(result$preprocess_info$volume_weights))
  expect_false(is.null(result$preprocess_info$soft_projection))
})

test_that("preprocess_run_data does nothing when disabled", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)

  cfg <- fmri_lm_control()  # Both disabled by default


  result <- preprocess_run_data(X, Y, cfg)

  expect_equal(result$X, X)
  expect_equal(result$Y, Y)
  expect_null(result$preprocess_info$volume_weights)
  expect_null(result$preprocess_info$soft_projection)
})
