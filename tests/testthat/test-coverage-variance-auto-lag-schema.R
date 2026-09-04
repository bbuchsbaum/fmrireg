# Fifteenth wave: HAC auto max-lag + result/variance schema helpers.

test_that(".fmri_lm_auto_max_lag covers run-list/null/censor and short-series edges", {
  skip_if_not_installed("fmriAR")
  set.seed(401)
  n <- 48L
  X <- cbind(1, rnorm(n), rnorm(n))
  runs <- rep(1:2, each = 24)

  lag <- fmrireg:::.fmri_lm_auto_max_lag(X, runs, cap = 12L)
  expect_true(is.integer(lag) || is.numeric(lag))
  expect_true(lag >= 0L && lag <= 12L)

  lag_list <- fmrireg:::.fmri_lm_auto_max_lag(
    X, runs = list(1:24, 25:48), censor = c(5L, 30L), cap = 8L
  )
  expect_true(lag_list >= 0L && lag_list <= 8L)

  lag_null <- fmrireg:::.fmri_lm_auto_max_lag(X, runs = NULL, cap = 6L)
  expect_true(lag_null >= 0L && lag_null <= 6L)

  # Short series forces zero lag (longest segment - 1 < 1 after cap)
  expect_equal(
    fmrireg:::.fmri_lm_auto_max_lag(X[1:2, , drop = FALSE], runs = NULL, cap = 60L),
    0L
  )
})

test_that(".fmri_lm_variance_from_context uses max_lag='auto' for HAC", {
  skip_if_not_installed("fmriAR")
  set.seed(402)
  n <- 40L
  V <- 2L
  X <- cbind(1, rnorm(n))
  E <- matrix(rnorm(n * V), n, V)
  ctx <- list(
    X = X,
    residuals = E,
    runs = rep(1:2, each = 20),
    censor = NULL,
    robust_weights = NULL
  )
  spec <- variance_spec(method = "hac", max_lag = "auto", taper = "tukey", df = "residual")
  expect_identical(spec$max_lag, "auto")

  fit <- suppressWarnings(fmrireg:::.fmri_lm_variance_from_context(ctx, spec))
  expect_equal(length(fit$covariance), V)
  expect_true(is.numeric(fit$max_lag))
  expect_true(fit$max_lag >= 0L)
  expect_equal(length(fit$df), V)
  expect_true(all(fit$df >= 1))
})

test_that("result schema helpers cover payload/df/weights/variance print", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  payload <- fmrireg:::.fmri_lm_beta_payload(fit$result)
  expect_true(is.null(payload) || inherits(payload, "tbl_df") || is.data.frame(payload))
  expect_null(fmrireg:::.fmri_lm_beta_payload(list()))
  expect_null(fmrireg:::.fmri_lm_beta_payload(list(betas = tibble::tibble())))

  nvox <- fmrireg:::.fmri_lm_result_nvox(fit$result)
  expect_true(is.numeric(nvox) && nvox >= 1)

  nominal <- fmrireg:::.fmri_lm_nominal_df(fit$result)
  expect_true(is.numeric(nominal))

  expect_equal(
    fmrireg:::.fmri_lm_df_from_weights(NULL, nominal = 20, nvox = 3L),
    rep(20, 3L)
  )
  dfs <- fmrireg:::.fmri_lm_df_from_weights(runif(30, 0.4, 1), nominal = 20, nvox = 2L)
  expect_equal(length(dfs), 2L)
  expect_true(all(dfs >= 1))

  # Per-voxel weight list path
  wlist <- list(runif(20, 0.5, 1), runif(20, 0.2, 0.8))
  dfs2 <- fmrireg:::.fmri_lm_df_from_weights(wlist, nominal = 18, nvox = 2L)
  expect_equal(length(dfs2), 2L)

  # Nested runwise weights fall back to nominal
  nested <- list(list(runif(10), runif(10)), list(runif(10), runif(10)))
  dfs3 <- fmrireg:::.fmri_lm_df_from_weights(nested, nominal = 15, nvox = 2L)
  expect_equal(dfs3, rep(15, 2L))

  vm <- fmrireg:::.new_fmri_variance_model(
    method = "model",
    covariance = diag(2),
    covariance_scope = "shared",
    scale = 1.2,
    df_nominal = 30,
    df_inference = 28,
    standard_errors = matrix(0.2, 3, 2)
  )
  expect_s3_class(vm, "fmri_lm_variance_model")
  expect_output(print(vm), "fmri_lm_variance_model|variance|model")
})
