# Closed-form validation of first-level variance.
#
# The package had extensive coverage of algebraic identities (F = t^2, F vs
# anova(), engine-vs-engine parity) but nothing that checked a reported standard
# error against an independently known covariance. These tests close that gap.
#
# Substrate (noise covariances, designs, the reduced-space oracle) lives in
# helper-variance.R and deliberately does not call fmrireg.

# ---------------------------------------------------------------- the oracle --
# Validate the reference implementation before using it to judge anything else.

test_that("reduced-space oracle reproduces L Sigma L' exactly", {
  n_t <- 200
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  P <- solve(crossprod(X), t(X))
  Cm <- rbind(c(1, 0, rep(0, p - 2)), c(1, -1, rep(0, p - 2)))
  L <- Cm %*% P
  run_id <- rep(1L, n_t)

  for (nm in names(.vt_regimes(n_t))) {
    Sigma <- .vt_regimes(n_t)[[nm]]
    got <- .vt_reduced_var(L, list(`1` = Sigma[1, ]), run_id)
    want <- L %*% Sigma %*% t(L)
    expect_equal(got, want, tolerance = 1e-10,
                 info = paste("covariance regime:", nm))
  }
})

test_that("oracle with Sigma = I reproduces the classical OLS formula", {
  n_t <- 160
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  Cm <- rbind(c(1, 0, rep(0, p - 2)), c(1, -1, rep(0, p - 2)))
  L <- Cm %*% solve(crossprod(X), t(X))

  got <- .vt_reduced_var(L, list(`1` = c(1, rep(0, 30))), rep(1L, n_t))
  want <- Cm %*% solve(crossprod(X)) %*% t(Cm)
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("oracle handles multi-run, run-specific, heteroscedastic noise", {
  n_run <- 3L; Tr <- 120L; n_t <- n_run * Tr
  run_id <- rep(seq_len(n_run), each = Tr)

  Xt <- do.call(rbind, lapply(seq_len(n_run),
                              function(r) .vt_design(Tr, "er")[, 1:2, drop = FALSE]))
  Nz <- matrix(0, n_t, 0)
  for (r in seq_len(n_run)) {
    D <- cbind(1, .vt_dct(Tr, 128))
    B <- matrix(0, n_t, ncol(D)); B[run_id == r, ] <- D
    Nz <- cbind(Nz, B)
  }
  X <- cbind(Xt, Nz); p <- ncol(X)
  L <- rbind(c(1, 0, rep(0, p - 2)), c(1, -1, rep(0, p - 2))) %*% solve(crossprod(X), t(X))

  # different noise structure AND different scale per run
  gl <- list(`1` = .vt_ar1(Tr, 0.5)[1, ] * 1.0,
             `2` = .vt_arma11(Tr, 0.6, 0.3)[1, ] * 1.4^2,
             `3` = .vt_1f(Tr)[1, ] * 0.8^2)
  Sigma <- matrix(0, n_t, n_t)
  for (r in seq_len(n_run)) Sigma[run_id == r, run_id == r] <- stats::toeplitz(gl[[r]])

  want <- L %*% Sigma %*% t(L)
  expect_equal(.vt_reduced_var(L, gl, run_id), want, tolerance = 1e-10)

  # Guard: ignoring run boundaries must be detected. Without this the test above
  # would pass for an implementation that concatenates runs.
  bad <- .vt_reduced_var(L, list(`1` = gl[[1]]), rep(1L, n_t))
  expect_gt(max(abs(bad - want)) / max(abs(want)), 1e-3)
})

# ------------------------------------------------- fmrireg vs closed form -----

test_that("beta_stats_matrix standard errors match summary(lm) exactly", {
  set.seed(101)
  n_t <- 180
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  y <- as.numeric(X %*% rnorm(p) + rnorm(n_t))

  fit <- lm(y ~ X - 1)
  se_lm <- unname(summary(fit)$coefficients[, "Std. Error"])

  XtXinv <- solve(crossprod(X))
  B <- matrix(coef(fit), ncol = 1)
  sigma <- sqrt(sum(residuals(fit)^2) / (n_t - p))
  res <- beta_stats_matrix(B, XtXinv, sigma, dfres = n_t - p,
                           varnames = paste0("v", seq_len(p)))

  expect_equal(as.numeric(res$data[[1]]$se[[1]]), se_lm, tolerance = 1e-10)
})

test_that(".fast_t_contrast standard error matches the closed form", {
  set.seed(102)
  n_t <- 180
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  y <- as.numeric(X %*% rnorm(p) + rnorm(n_t))
  fit <- lm(y ~ X - 1)

  XtXinv <- solve(crossprod(X))
  B <- matrix(coef(fit), ncol = 1)
  s2 <- sum(residuals(fit)^2) / (n_t - p)
  cvec <- c(1, -1, rep(0, p - 2))

  out <- .fast_t_contrast(B, sigma2 = s2, XtXinv = XtXinv,
                          l = matrix(cvec, nrow = 1), df = n_t - p)

  expect_equal(as.numeric(out$se), .vt_iid_se(X, s2, cvec), tolerance = 1e-10)
  # and against an independent route through lm()
  expect_equal(as.numeric(out$se),
               unname(sqrt(t(cvec) %*% vcov(fit) %*% cvec)[1, 1]),
               tolerance = 1e-10)
})

test_that("fmri_lm reports OLS standard errors that match lm()", {
  skip_if_not_installed("fmridataset")
  set.seed(103)
  # single run: runwise fitting then coincides with one concatenated fit, so the
  # comparison isolates the SE arithmetic from the run-pooling scheme
  TR <- 2; Tr <- 120L; n_run <- 1L; n_t <- Tr * n_run

  onsets <- sort(runif(30, 0, Tr * TR - 20))
  cond <- factor(rep(c("a", "b"), length.out = 30))
  etab <- data.frame(onset = rep(onsets, n_run),
                     cond = rep(cond, n_run),
                     run = rep(seq_len(n_run), each = 30))
  Y <- matrix(rnorm(n_t * 3), n_t, 3)
  dset <- matrix_dataset(Y, TR = TR, run_length = rep(Tr, n_run), event_table = etab)

  fit <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                 strategy = "runwise", ar_options = list(struct = "iid"))
  se_pkg <- as.matrix(standard_error(fit, type = "estimates"))

  # Independent route: recover the design fmrireg built and refit with lm()
  X <- as.matrix(design_matrix(fit$model))
  se_lm <- vapply(seq_len(ncol(Y)), function(v) {
    unname(summary(lm(Y[, v] ~ X - 1))$coefficients[1, "Std. Error"])
  }, numeric(1))

  expect_equal(unname(se_pkg[, 1]), se_lm, tolerance = 1e-6)
})

# ------------------------------------------- external reference agreement -----

test_that("calculate_sandwich_variance matches sandwich::vcovHC(HC0)", {
  skip_if_not_installed("sandwich")
  set.seed(104)
  n_t <- 150
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  # heteroskedastic errors, so the sandwich genuinely differs from the OLS vcov
  y <- as.numeric(X %*% rnorm(p) + rnorm(n_t) * seq(0.5, 2, length.out = n_t))
  fit <- lm(y ~ X - 1)
  r <- matrix(residuals(fit), ncol = 1)

  # the function applies an n/(n-p) small-sample correction, i.e. HC1
  got <- calculate_sandwich_variance(X, r, solve(crossprod(X)))
  want <- sandwich::vcovHC(fit, type = "HC1")

  expect_equal(unname(sqrt(diag(got))), unname(sqrt(diag(want))), tolerance = 1e-8)
})

test_that("compute_sandwich_variance matches sandwich::vcovHC", {
  skip_if_not_installed("sandwich")
  set.seed(105)
  n_t <- 150
  X <- .vt_design(n_t, "er")
  p <- ncol(X)
  y <- as.numeric(X %*% rnorm(p) + rnorm(n_t) * seq(0.5, 2, length.out = n_t))
  fit <- lm(y ~ X - 1)
  r <- matrix(residuals(fit), ncol = 1)

  for (ty in c("HC0", "HC1")) {
    got <- compute_sandwich_variance(X, r, type = ty)
    want <- sandwich::vcovHC(fit, type = ty)
    expect_equal(unname(sqrt(diag(got))), unname(sqrt(diag(want))),
                 tolerance = 1e-8, info = ty)
  }
})

test_that("bias correction recovers the true autocovariance from residuals", {
  n_t <- 220
  X <- .vt_design(n_t, "er")
  H <- 25L
  A <- .vt_bias_map(X, H)
  M <- diag(n_t) - X %*% solve(crossprod(X), t(X))
  d <- abs(outer(seq_len(n_t), seq_len(n_t), "-"))

  # AR and ARMA covariances are supported within the lag budget; 1/f and the
  # oscillatory regime are not, and are excluded by design (see below).
  for (nm in c("AR(1) rho=0.4", "ARMA(1,1)")) {
    Sigma <- .vt_regimes(n_t)[[nm]]
    g_true <- Sigma[1, seq_len(H + 1L)]
    g_raw <- .vt_lag_mean(M %*% Sigma %*% M, H, d)
    g_corr <- solve(A, g_raw)

    err_raw <- max(abs(g_raw / g_raw[1] - g_true / g_true[1]))
    err_cor <- max(abs(g_corr / g_corr[1] - g_true / g_true[1]))

    # the correction must work ...
    expect_lt(err_cor, 1e-6)
    # ... and the uncorrected estimator must be measurably biased, otherwise
    # the assertion above is vacuous and would pass for a no-op correction.
    expect_gt(err_raw, 0.02)
  }
})

test_that("AR effective df is no longer phi-blind or deflated", {
  # This began as a change detector pinning the old behaviour: calculate_effective_df()
  # multiplied n - p by 1/(1 + 2*sum(1 - k/n)) whenever ar_order > 0, returning 62.9
  # here regardless of the true autocorrelation, on top of whitening. Fixed
  # 2026-08-07; the detector has done its job and now asserts the corrected
  # contract. Full coverage lives in test-effective-df-correctness.R.
  n <- 200; p <- 12

  expect_equal(calculate_effective_df(n, p, ar_order = 1), n - p)
  expect_equal(calculate_effective_df(n, p, ar_order = 2), n - p)

  # the specific deflated value must not come back
  expect_false(isTRUE(all.equal(calculate_effective_df(n, p, ar_order = 1),
                                62.9, tolerance = 1e-2)))

  # "satterthwaite" is a real method now, not an alias for the simple one
  expect_error(calculate_effective_df(n, p, ar_order = 1, method = "satterthwaite"),
               "requires both")
})

test_that("residual bias correction is limited by the lag budget", {
  # Documents a real boundary of the method rather than asserting it away: when
  # the true covariance has support beyond lag H, the truncated correction
  # system cannot recover it. Consumers must either widen H or high-pass filter.
  n_t <- 220
  X <- .vt_design(n_t, "er")
  H <- 25L
  A <- .vt_bias_map(X, H)
  M <- diag(n_t) - X %*% solve(crossprod(X), t(X))
  d <- abs(outer(seq_len(n_t), seq_len(n_t), "-"))

  Sigma <- .vt_regimes(n_t)[["AR+oscillatory"]]
  g_true <- Sigma[1, seq_len(H + 1L)]
  g_corr <- solve(A, .vt_lag_mean(M %*% Sigma %*% M, H, d))
  err <- max(abs(g_corr / g_corr[1] - g_true / g_true[1]))

  expect_gt(err, 1e-6)   # known limitation, asserted so a future fix is noticed
  expect_lt(err, 0.2)    # but it must not be catastrophic
})
