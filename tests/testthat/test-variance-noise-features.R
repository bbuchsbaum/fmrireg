# Integration between fmrireg's design-side variance machinery and fmriAR's
# noise-covariance features (fmriAR >= 0.3.3: noise_acvf(), acvf_bias_matrix()).
#
# The split under test:
#   fmriAR  owns Sigma  -- autocovariance, run segmentation, censoring, debiasing
#   fmrireg owns L      -- the design/contrast projection and the combination
# so Cov(C beta) = sum_h gamma_h (L T_h L').
#
# Oracle (helper-variance.R) calls neither package, so agreement between the two
# is genuine cross-implementation evidence rather than a tautology.

skip_if_no_acvf <- function() {
  skip_if_not_installed("fmriAR")
  if (!all(c("noise_acvf", "acvf_bias_matrix") %in% getNamespaceExports("fmriAR"))) {
    skip("fmriAR >= 0.3.3 required for noise_acvf()/acvf_bias_matrix()")
  }
}

test_that("fmriAR's residual-bias matrix matches the independent oracle", {
  skip_if_no_acvf()
  n_t <- 240L
  H <- 20L
  for (ty in c("er", "block")) {
    X <- .vt_design(n_t, ty)
    got <- fmriAR::acvf_bias_matrix(X, max_lag = H)[[1]]
    want <- .vt_bias_map(X, H)
    expect_equal(unname(got), unname(want), tolerance = 1e-8, info = ty)
  }
})

test_that("lag budget is governed by the conditioning of the bias matrix", {
  # cond(A) depends only on the design, so a safe max_lag can be chosen before
  # any data is seen. It grows steeply, and past roughly 1e2 the debiasing solve
  # amplifies noise faster than the extra lags buy accuracy.
  skip_if_no_acvf()
  n_t <- 300L
  X <- .vt_design(n_t, "er")

  kap <- vapply(c(15L, 25L, 40L, 60L), function(H) {
    kappa(fmriAR::acvf_bias_matrix(X, max_lag = H)[[1]], exact = TRUE)
  }, numeric(1))

  expect_true(all(diff(kap) > 0))       # monotone in max_lag
  expect_lt(kap[1], 50)                 # short budgets are well conditioned
  expect_gt(kap[4], 1e3)                # long budgets are not
})

test_that("noise_acvf recovers autocovariance shape and scale after debiasing", {
  skip_if_no_acvf()
  set.seed(201)
  n_t <- 300L
  H <- 20L
  V <- 2000L
  X <- .vt_design(n_t, "er")
  M <- diag(n_t) - X %*% solve(crossprod(X), t(X))
  scale2 <- 2.5^2

  for (nm in c("AR(1) rho=0.4", "ARMA(1,1)")) {
    Sigma <- .vt_regimes(n_t)[[nm]] * scale2
    E <- crossprod(chol(Sigma), matrix(rnorm(n_t * V), n_t, V))
    R <- M %*% E
    g_true <- Sigma[1, seq_len(H + 1L)]

    raw <- fmriAR::noise_acvf(R, max_lag = H)
    cor <- fmriAR::noise_acvf(R, max_lag = H, design = X, correction_max_lag = H)
    expect_false(isTRUE(raw$corrected))
    expect_true(isTRUE(cor$corrected))

    g_raw <- as.numeric(raw$acvf[[1]])
    g_cor <- as.numeric(cor$acvf[[1]])
    err <- function(g) max(abs(g / g[1] - g_true / g_true[1]))

    expect_lt(err(g_cor), 0.02)                       # shape recovered
    expect_lt(abs(g_cor[1] / g_true[1] - 1), 0.02)    # and the scale, not just the shape
    expect_lt(err(g_cor), err(g_raw))                 # correction is not a no-op
    expect_gt(err(g_raw), 0.05)                       # ... and had something to fix
  }
})

test_that("noise_acvf segments at run boundaries and censored frames", {
  skip_if_no_acvf()
  set.seed(202)
  Tr <- 150L; n_run <- 2L; n_t <- Tr * n_run
  runs <- rep(seq_len(n_run), each = Tr)
  V <- 1500L

  # white within run, large level shift between runs: a run-blind estimator
  # reads the shift as long-range autocorrelation
  E <- matrix(rnorm(n_t * V), n_t, V)
  E[runs == 2, ] <- E[runs == 2, ] + 6

  blind <- as.numeric(fmriAR::noise_acvf(E, max_lag = 6L)$acvf[[1]])
  aware <- as.numeric(fmriAR::noise_acvf(E, runs = runs, max_lag = 6L)$acvf[[1]])

  expect_gt(abs(blind[2] / blind[1]), 0.5)     # the bug this guards against
  expect_lt(abs(aware[2] / aware[1]), 0.05)    # correct behaviour

  # censored frames must not contribute
  E2 <- matrix(rnorm(n_t * V), n_t, V)
  bad <- c(40:45, 200:205)
  E2[bad, ] <- E2[bad, ] + 50
  keep <- as.numeric(fmriAR::noise_acvf(E2, runs = runs, max_lag = 6L)$acvf[[1]])
  drop <- as.numeric(fmriAR::noise_acvf(E2, runs = runs, censor = bad,
                                        max_lag = 6L)$acvf[[1]])
  expect_gt(keep[1], 10)                       # corrupted frames dominate the variance
  expect_lt(abs(drop[1] - 1), 0.1)             # ... and are excluded when censored
})

test_that("end-to-end: fmriAR features + reduced-space kernels give correct SEs", {
  # The whole point of the split. Estimate gamma with fmriAR, project it through
  # the design with the reduced-space identity, and check the resulting contrast
  # standard error against the analytically known sqrt(L Sigma L').
  skip_if_no_acvf()
  set.seed(203)
  n_t <- 300L
  H <- 20L
  V <- 400L
  run_id <- rep(1L, n_t)

  for (ty in c("er", "block")) {
    X <- .vt_design(n_t, ty)
    p <- ncol(X)
    P <- solve(crossprod(X), t(X))
    cvec <- c(1, -1, rep(0, p - 2))
    L <- matrix(cvec, nrow = 1) %*% P

    for (nm in c("AR(1) rho=0.4", "ARMA(1,1)")) {
      Sigma <- .vt_regimes(n_t)[[nm]] * 1.7^2
      Y <- crossprod(chol(Sigma), matrix(rnorm(n_t * V), n_t, V))
      R <- Y - X %*% (P %*% Y)

      g <- fmriAR::noise_acvf(R, max_lag = H, design = X,
                              correction_max_lag = H)$acvf[[1]]
      se_hat <- sqrt(drop(.vt_reduced_var(L, list(`1` = as.numeric(g)), run_id)))
      se_true <- .vt_true_se(X, Sigma, cvec)

      expect_lt(abs(se_hat / se_true - 1), 0.05,
                label = sprintf("%s / %s: relative SE error", ty, nm))

      # the naive iid standard error must be materially wrong here, otherwise
      # the check above would pass for an estimator that ignores the noise model
      s2_naive <- sum(R[, 1]^2) / (n_t - p)
      se_naive <- .vt_iid_se(X, s2_naive, cvec)
      expect_gt(abs(se_naive / se_true - 1), 0.05)
    }
  }
})
