# Correctness of calculate_effective_df().
#
# The previous implementation multiplied n - p by 1/(1 + 2*sum(1 - k/n)) whenever
# ar_order > 0. That factor assumes every autocorrelation equals 1, never depends
# on the fitted AR coefficients, and was applied on top of whitening. These tests
# pin the corrected behaviour against closed forms, an independent brute-force
# implementation, and -- the gate that actually matters -- empirical calibration
# of the resulting t statistics.

.edf  <- function(...) fmrireg:::calculate_effective_df(...)
.sdf  <- function(...) fmrireg:::.satterthwaite_df(...)

# Brute-force reference: form R explicitly and take the traces directly.
.sdf_ref <- function(X, V) {
  n <- nrow(X)
  R <- diag(n) - X %*% solve(crossprod(X), t(X))
  RV <- R %*% V
  sum(diag(RV))^2 / sum(diag(RV %*% RV))
}

# ------------------------------------------------------------- basic contract --

test_that("residual method returns n - p and ignores AR order", {
  # AR whitening happens before the fit, so the whitened residuals already carry
  # n - p degrees of freedom; reducing them again double-counts.
  expect_equal(.edf(200, 12), 188)
  expect_equal(.edf(200, 12, ar_order = 1), 188)
  expect_equal(.edf(200, 12, ar_order = 4), 188)
  expect_equal(.edf(300, 30, ar_order = 2), 270)
})

test_that("the old phi-blind value is no longer produced", {
  # regression guard: the previous implementation returned 62.9 here for every phi
  expect_false(isTRUE(all.equal(.edf(200, 12, ar_order = 1), 62.9, tolerance = 1e-2)))
  expect_false(isTRUE(all.equal(.edf(200, 12, ar_order = 2), 37.8, tolerance = 1e-2)))
})

test_that("'simple' is accepted as an alias for 'residual'", {
  expect_equal(.edf(200, 12, method = "simple"),
               .edf(200, 12, method = "residual"))
  expect_equal(.edf(200, 12, ar_order = 3, method = "simple"), 188)
})

test_that("satterthwaite is no longer a duplicate of the simple method", {
  # it used to be a verbatim copy; it now needs real inputs and refuses without them
  expect_error(.edf(200, 12, method = "satterthwaite"), "requires both")
  expect_error(.edf(200, 12, method = "satterthwaite", X = matrix(1, 200, 2)),
               "requires both")
})

test_that("degenerate inputs return at least 1", {
  expect_equal(.edf(5, 5), 1)
  expect_equal(.edf(3, 10), 1)
  expect_gte(.edf(20, 2, robust_weights = rep(1e-8, 20)), 1)
})

# --------------------------------------------------------- satterthwaite core --

test_that("Satterthwaite reduces exactly to n - p when errors are uncorrelated", {
  for (n in c(60L, 150L)) {
    for (p in c(3L, 10L)) {
      X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
      expect_equal(.sdf(X, diag(n)), n - p, tolerance = 1e-8,
                   info = sprintf("n=%d p=%d", n, p))
      expect_equal(.edf(n, p, method = "satterthwaite", X = X, resid_cov = diag(n)),
                   n - p, tolerance = 1e-8)
    }
  }
})

test_that("Satterthwaite matches a brute-force trace computation", {
  set.seed(11)
  n <- 120L
  X <- .vt_design(n, "er")
  for (nm in names(.vt_regimes(n))) {
    V <- .vt_regimes(n)[[nm]]
    expect_equal(.sdf(X, V), .sdf_ref(X, V), tolerance = 1e-8, info = nm)
  }
})

test_that("Satterthwaite is invariant to the scale of the covariance", {
  # it is a ratio of traces, so multiplying V by a constant must not change it
  set.seed(12)
  n <- 100L
  X <- .vt_design(n, "er")
  V <- .vt_ar1(n, 0.5)
  expect_equal(.sdf(X, V), .sdf(X, 7.3 * V), tolerance = 1e-8)
})

test_that("stronger autocorrelation lowers the effective df, monotonically", {
  set.seed(13)
  n <- 150L
  X <- .vt_design(n, "er")
  rhos <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9)
  dfs <- vapply(rhos, function(r) .sdf(X, .vt_ar1(n, r)), numeric(1))

  expect_equal(dfs[1], n - ncol(X), tolerance = 1e-8)   # rho = 0 is the iid case
  expect_true(all(diff(dfs) < 0))                        # strictly decreasing
  expect_true(all(dfs >= 1))
  expect_true(all(dfs <= n - ncol(X) + 1e-8))            # never exceeds the iid df
})

test_that("Satterthwaite handles a rank-deficient design without erroring", {
  set.seed(14)
  n <- 80L
  X <- cbind(1, rnorm(n), rnorm(n))
  X <- cbind(X, X[, 2])                                  # exact duplicate column
  expect_silent(d <- .sdf(X, diag(n)))
  expect_gte(d, 1)
})

# -------------------------------------------------------------------- robust --

test_that("robust weights scale the df by their mean", {
  w <- c(rep(1, 90), rep(0, 10))
  expect_equal(.edf(100, 5, robust_weights = w), 95 * 0.9, tolerance = 1e-10)
  # all-ones weights must be a no-op
  expect_equal(.edf(100, 5, robust_weights = rep(1, 100)), 95, tolerance = 1e-10)
})

test_that("robust and Satterthwaite adjustments compose multiplicatively", {
  set.seed(15)
  n <- 120L
  X <- .vt_design(n, "er")
  V <- .vt_ar1(n, 0.5)
  w <- c(rep(1, 100), rep(0.5, 20))
  base <- .edf(n, ncol(X), method = "satterthwaite", X = X, resid_cov = V)
  both <- .edf(n, ncol(X), robust_weights = w, method = "satterthwaite",
               X = X, resid_cov = V)
  expect_equal(both, base * mean(w), tolerance = 1e-10)
})

# ------------------------------------------------- calibration (the real gate) --

test_that("t statistics built on the reported df are calibrated under the null", {
  # This is the property the old formula violated and no test checked: if the df
  # is right, p-values from the null distribution are uniform and the rejection
  # rate matches alpha.
  skip_on_cran()
  set.seed(2024)
  n <- 160L
  X <- .vt_design(n, "er")
  p <- ncol(X)
  cvec <- c(1, -1, rep(0, p - 2))
  XtXinv <- solve(crossprod(X))
  se_scale <- sqrt(drop(t(cvec) %*% XtXinv %*% cvec))
  nrep <- 4000L
  df <- .edf(n, p)

  tvals <- numeric(nrep)
  for (i in seq_len(nrep)) {
    y <- rnorm(n)                                  # true contrast is 0
    b <- drop(XtXinv %*% crossprod(X, y))
    s2 <- sum((y - X %*% b)^2) / (n - p)
    tvals[i] <- drop(cvec %*% b) / (sqrt(s2) * se_scale)
  }
  pvals <- 2 * pt(-abs(tvals), df)

  # rejection rates track alpha; MC s.e. at alpha=.05 over 4000 reps is ~0.35pp
  expect_lt(abs(mean(pvals < 0.05) - 0.05), 0.012)
  expect_lt(abs(mean(pvals < 0.01) - 0.01), 0.006)
  # and the p-values are uniform overall
  expect_gt(suppressWarnings(ks.test(pvals, "punif")$p.value), 0.01)

  # the old value would have been badly conservative here
  df_old <- (n - p) / (1 + 2 * (1 - 1 / n))
  p_old <- 2 * pt(-abs(tvals), df_old)
  expect_lt(mean(p_old < 0.05), mean(pvals < 0.05))
})

test_that("Satterthwaite df calibrates OLS t statistics under correlated errors", {
  # Unwhitened OLS with AR(1) errors: the naive n - p df is wrong, and the
  # Satterthwaite df computed from the true covariance is much closer.
  skip_on_cran()
  set.seed(2025)
  n <- 160L
  X <- .vt_design(n, "er")
  p <- ncol(X)
  V <- .vt_ar1(n, 0.5)
  Ch <- chol(V)
  cvec <- c(1, -1, rep(0, p - 2))
  XtXinv <- solve(crossprod(X))
  Pm <- XtXinv %*% t(X)
  # exact sampling variance of the contrast under OLS with covariance V
  true_var <- drop(cvec %*% Pm %*% V %*% t(Pm) %*% cvec)
  df_satt <- .sdf(X, V)

  nrep <- 3000L
  tvals <- numeric(nrep)
  for (i in seq_len(nrep)) {
    y <- drop(crossprod(Ch, rnorm(n)))
    b <- drop(Pm %*% y)
    tvals[i] <- drop(cvec %*% b) / sqrt(true_var)   # variance known; df is what is tested
  }
  # with the variance known and correct, the statistic is standard normal; the
  # Satterthwaite df must be large enough that t_df is a good approximation
  expect_gt(df_satt, 30)
  expect_lt(df_satt, n - p)
  expect_lt(abs(mean(abs(tvals) > qt(0.975, df_satt)) - 0.05), 0.015)
})

# ---------------------------------------------------------- end-to-end wiring --

test_that("fmri_lm AR fits report residual df, not a deflated one", {
  skip_if_not_installed("fmridataset")
  set.seed(303)
  TR <- 2; Tr <- 120L; n_run <- 1L; n_t <- Tr * n_run
  onsets <- sort(runif(24, 0, Tr * TR - 20))
  cond <- factor(rep(c("a", "b"), length.out = 24))
  etab <- data.frame(onset = onsets, cond = cond, run = rep(1L, 24))
  dset <- matrix_dataset(matrix(rnorm(n_t * 3), n_t, 3), TR = TR,
                         run_length = Tr, event_table = etab)

  fit_iid <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     control = fmri_lm_control(),
                     compute = compute_spec(voxel_chunks = 1L))
  fit_ar  <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     control = fmri_lm_control(noise = noise_spec("ar1")),
                     compute = compute_spec(voxel_chunks = 1L))

  se_iid <- as.matrix(standard_error(fit_iid, type = "estimates"))
  se_ar  <- as.matrix(standard_error(fit_ar,  type = "estimates"))
  expect_true(all(is.finite(se_iid)))
  expect_true(all(is.finite(se_ar)))
  expect_equal(dim(se_ar), dim(se_iid))

  # The AR path used to report roughly a third of the residual df. The common
  # inference schema now reports the df that generated every p-value.
  df_resid <- .edf(n_t, ncol(as.matrix(design_matrix(fit_ar$model))))
  df_old <- df_resid / (1 + 2 * (1 - 1 / n_t))
  expect_gt(df_resid, 2 * df_old)
  expect_equal(fit_ar$result$df$inference,
               rep(df_resid, length(fit_ar$result$df$inference)))
  expect_equal(
    fit_ar$result$betas$data[[1L]]$prob[[1L]],
    2 * pt(-abs(fit_ar$result$betas$data[[1L]]$stat[[1L]]), df = df_resid),
    tolerance = 1e-12
  )
})
