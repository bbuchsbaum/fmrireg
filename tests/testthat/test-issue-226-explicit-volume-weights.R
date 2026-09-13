make_issue_226_fixture <- function(nvox = 1L) {
  set.seed(2)
  n <- 80L
  events <- data.frame(
    onset = seq(6, 132, by = 18),
    condition = factor(rep(c("A", "B"), 4)),
    run = 1L
  )
  sframe <- sampling_frame(n, TR = 2)
  event_design <- as.matrix(design_matrix(event_model(
    onset ~ hrf(condition), data = events, block = ~run,
    sampling_frame = sframe
  )))
  # fmri_model() orders event columns before baseline columns.
  full_design <- cbind(event_design, `(Intercept)` = 1)
  event_columns <- seq_len(ncol(event_design))
  signal <- drop(full_design %*% c(2, -1, 100))
  Y <- vapply(seq_len(nvox), function(v) {
    signal + rnorm(n, sd = 0.2)
  }, numeric(n))
  bad <- c(10L, 11L, 12L, 40L)
  Y[bad, ] <- Y[bad, , drop = FALSE] + 50
  weights <- as.numeric(!(seq_len(n) %in% bad))
  dataset <- matrix_frame(
    Y, TR = 2, run_length = n, event_table = events
  )
  baseline <- baseline_model(basis = "constant", sframe = sframe)

  list(
    Y = Y,
    full_design = full_design,
    event_columns = event_columns,
    weights = weights,
    dataset = dataset,
    baseline = baseline
  )
}

fit_issue_226 <- function(fixture, weights, chunks = 1L,
                          backend = "matrix") {
  fmri_lm(
    onset ~ hrf(condition), block = ~run,
    dataset = fixture$dataset,
    baseline_model = fixture$baseline,
    control = fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec("iid"),
      weights = weights
    ),
    compute = compute_spec(voxel_chunks = chunks, backend = backend)
  )
}

test_that("explicit values enable volume weighting", {
  values <- c(1, 0, 0.5)
  spec <- weights_spec(values = values)

  expect_true(spec$enabled)
  expect_identical(spec$method, "none")
  expect_equal(spec$values, values)
  expect_equal(spec$weights, values)
  expect_error(weights_spec(values = c(1, -1)), "nonnegative")
})

test_that("joint IID OLS honors explicit zero weights", {
  fixture <- make_issue_226_fixture()
  keep <- fixture$weights > 0
  oracle <- lm.fit(
    fixture$full_design[keep, , drop = FALSE],
    fixture$Y[keep, , drop = FALSE]
  )
  oracle_coefficients <- as.matrix(oracle$coefficients)

  fit <- expect_no_warning(fit_issue_226(
    fixture, weights_spec(values = fixture$weights)
  ))

  expect_equal(
    as.numeric(coef(fit)),
    as.numeric(oracle_coefficients[fixture$event_columns, , drop = FALSE]),
    tolerance = 1e-10
  )
  expect_equal(fit$result$rdf, oracle$df.residual)
  expect_equal(
    as.numeric(fit$result$covariance_model_basis),
    as.numeric(solve(crossprod(fixture$full_design[keep, , drop = FALSE]))),
    tolerance = 1e-10
  )
  expected_rss <- sum(oracle$residuals^2)
  expect_equal(fit$result$resvar, expected_rss / oracle$df.residual,
               tolerance = 1e-10)
  expect_equal(
    fit$result$cov.unscaled,
    fit$result$covariance_model_basis,
    tolerance = 0
  )
})

test_that("joint IID weighted estimates are invariant to voxel chunking", {
  fixture <- make_issue_226_fixture(nvox = 3L)
  fractional <- seq(0.25, 1.75, length.out = nrow(fixture$Y))
  oracle <- lm.wfit(fixture$full_design, fixture$Y, w = fractional)

  one <- fit_issue_226(fixture, weights_spec(values = fractional), chunks = 1L)
  three <- fit_issue_226(
    fixture,
    weights_spec("inverse_squared", values = fractional),
    chunks = 3L
  )

  expected <- as.matrix(oracle$coefficients)[fixture$event_columns, , drop = FALSE]
  expect_equal(unname(as.matrix(coef(one))), unname(expected), tolerance = 1e-10)
  expect_equal(unname(as.matrix(coef(three))), unname(expected), tolerance = 1e-10)
  expect_equal(coef(three), coef(one), tolerance = 1e-12)
  expect_equal(one$result$rdf, oracle$df.residual)
})

test_that("unsupported weighted backends fail instead of fitting unweighted", {
  fixture <- make_issue_226_fixture()

  expect_error(
    fit_issue_226(fixture, weights_spec("tukey")),
    "Computed volume weights are not supported by joint IID OLS"
  )
  expect_error(
    fit_issue_226(
      fixture, weights_spec(values = fixture$weights), backend = "reference"
    ),
    "Volume weighting is not supported by the reference fitting backend"
  )
  expect_error(
    fmri_lm(
      onset ~ hrf(condition), block = ~run,
      dataset = fixture$dataset,
      baseline_model = fixture$baseline,
      control = fmri_lm_control(
        estimation = estimation_spec("joint"),
        noise = noise_spec("ar1"),
        weights = weights_spec(values = fixture$weights)
      ),
      compute = compute_spec(voxel_chunks = 1L)
    ),
    "Volume weighting is not supported by joint AR or robust fitting"
  )
})
