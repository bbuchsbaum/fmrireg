# Review 2026-10, stream D, item P-D3: HAC / sandwich / Satterthwaite inference.
# The per-voxel work is checked against sandwich::vcovHAC on lm fits, and the
# design-only work (bread) is shown to be done once.

.vd_fixture <- function(n = 120L, V = 4L, seed = 3L) {
  set.seed(seed)
  ev <- data.frame(onset = seq(6, 2 * n - 30, by = 14),
                   condition = factor(rep(c("A", "B"), length.out = length(seq(6, 2 * n - 30, by = 14)))),
                   run = 1)
  Y <- vapply(seq_len(V), function(v) as.numeric(arima.sim(list(ar = 0.4), n)) + 20,
              numeric(n))
  list(dset = matrix_frame(Y, TR = 2, run_length = n, event_table = ev), Y = Y)
}

.vd_fit <- function(fx, variance) {
  fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$dset, durations = 0,
          control = fmri_lm_control(estimation = estimation_spec("joint"),
                                    variance = variance),
          compute = compute_spec(voxel_chunks = 2))
}

test_that("score covariance equals sandwich::vcovHAC for each taper", {
  skip_if_not_installed("sandwich")
  set.seed(9)
  n <- 150
  X <- cbind(1, rnorm(n), sin(seq_len(n) / 7))
  for (v in 1:3) {
    y <- drop(X %*% c(1, 0.5, -1)) + as.numeric(arima.sim(list(ar = 0.5), n))
    fit <- lm(y ~ X - 1)
    for (taper in c("none", "tukey", "parzen")) {
      L <- 4L
      w <- vapply(0:L, function(h) fmrireg:::.fmri_lm_hac_weight(h, L, taper),
                  numeric(1))
      ref <- sandwich::vcovHAC(fit, weights = w, prewhite = FALSE, adjust = TRUE)
      got <- fmrireg:::.fmri_lm_score_covariance(
        X, residuals(fit), list(seq_len(n)), max_lag = L, taper = taper,
        debias = TRUE
      )
      expect_equal(unname(got), unname(ref), tolerance = 1e-8,
                   info = paste(taper, v))
    }
    # max_lag = 0 is the HC0 sandwich with the n / (n - p) correction (HC1).
    got0 <- fmrireg:::.fmri_lm_score_covariance(X, residuals(fit), list(seq_len(n)))
    expect_equal(unname(got0), unname(sandwich::vcovHC(fit, type = "HC1")),
                 tolerance = 1e-8)
  }
})

test_that("fmri_lm HAC standard errors match vcovHAC voxel by voxel", {
  skip_if_not_installed("sandwich")
  fx <- .vd_fixture()
  fit <- .vd_fit(fx, variance_spec("hac", max_lag = 3, taper = "parzen"))
  X <- as.matrix(design_matrix(fit$model))
  se <- fit$result$betas$data[[1]]$se[[1]]
  w <- vapply(0:3, function(h) fmrireg:::.fmri_lm_hac_weight(h, 3L, "parzen"),
              numeric(1))
  for (v in seq_len(ncol(fx$Y))) {
    ref <- sandwich::vcovHAC(lm(fx$Y[, v] ~ X - 1), weights = w,
                             prewhite = FALSE, adjust = TRUE)
    expect_equal(unname(se[v, ]), unname(sqrt(diag(ref))), tolerance = 1e-8)
  }
})

test_that("model-based Satterthwaite inference builds no sandwich covariance", {
  fx <- .vd_fixture(V = 6L)
  score_calls <- 0L
  orig_score <- fmrireg:::.fmri_lm_score_covariance
  testthat::local_mocked_bindings(
    .fmri_lm_score_covariance = function(...) {
      score_calls <<- score_calls + 1L
      orig_score(...)
    }
  )
  m <- .vd_fit(fx, variance_spec("model", df = "satterthwaite"))
  expect_equal(score_calls, 0L)

  # HAC needs exactly one covariance per voxel.
  h <- .vd_fit(fx, variance_spec("hac", max_lag = 2))
  expect_equal(score_calls, ncol(fx$Y))
})

test_that("the sandwich bread is computed once per fit, not per voxel", {
  fx <- .vd_fixture(V = 6L)
  bread_calls <- 0L
  orig_bread <- fmrireg:::.fmri_lm_bread
  testthat::local_mocked_bindings(
    .fmri_lm_bread = function(...) {
      bread_calls <<- bread_calls + 1L
      orig_bread(...)
    }
  )
  h <- .vd_fit(fx, variance_spec("hac", max_lag = 2))
  expect_equal(bread_calls, 1L)
})
