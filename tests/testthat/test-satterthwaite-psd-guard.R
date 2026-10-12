# Performance guard for the PSD-projected Satterthwaite df: voxels whose
# covariance the projection changes are handled by the banded eigensolver
# (src/satterthwaite_psd.cpp), which forms no n x n or segment-sized square
# matrix and calls no eigen(). The dense per-block helper is only the fallback
# for a failed LAPACK call.

.spg_problem <- function(n = 240L, p = 6L, seed = 2L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  E <- vapply(c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95), function(a) {
    as.numeric(arima.sim(list(ar = a), n))
  }, numeric(n))
  E <- E - X %*% qr.solve(X, E)
  list(X = X, E = E, runs = list(1:120, 121:240))
}

.spg_df <- function(pr, max_lag = 30L, taper = "none") {
  ctx <- list(X = pr$X, residuals = pr$E, runs = pr$runs, censor = NULL)
  fmrireg:::.fmri_lm_variance_from_context(
    ctx, variance_spec("model", max_lag = max_lag, taper = taper,
                       df = "satterthwaite")
  )$df
}

.spg_n_projected <- function(pr, max_lag = 30L, taper = "none") {
  segs <- fmrireg:::.fmri_lm_segments(nrow(pr$X), pr$runs)
  gl <- fmrireg:::.fmri_lm_segment_gammas(pr$E, segs, max_lag, taper)
  sum(fmrireg:::.fmri_lm_toeplitz_needs_projection(gl$gamma, gl$lags,
                                                   gl$offsets))
}

test_that("projected Satterthwaite df calls no eigen()", {
  pr <- .spg_problem()
  # Every voxel here needs the exact projection (truncated window, lag 30).
  expect_equal(.spg_n_projected(pr), ncol(pr$E))
  n_eigen <- 0L
  suppressMessages(trace("eigen", where = baseenv(), print = FALSE,
                         tracer = function() n_eigen <<- n_eigen + 1L))
  on.exit(suppressMessages(untrace("eigen", where = baseenv())), add = TRUE)
  df <- .spg_df(pr)
  expect_length(df, ncol(pr$E))
  expect_true(all(is.finite(df)))
  expect_equal(n_eigen, 0L)
})

test_that("projected path is one banded call with no dense fallback", {
  pr <- .spg_problem()
  calls <- 0L
  statuses <- integer(0)
  clipped <- integer(0)
  orig_cpp <- fmrireg:::satterthwaite_psd_traces_cpp
  testthat::local_mocked_bindings(
    satterthwaite_psd_traces_cpp = function(...) {
      calls <<- calls + 1L
      res <- orig_cpp(...)
      statuses <<- c(statuses, res$status)
      clipped <<- c(clipped, res$n_clipped)
      res
    },
    .fmri_lm_projected_traces_dense = function(...) {
      stop("dense per-block fallback used")
    },
    .fmri_lm_psd = function(...) stop("dense PSD projection used"),
    .fmri_lm_residual_covariance = function(...) {
      stop("dense residual covariance formed")
    }
  )
  for (cfg in list(list(30L, "none"), list(8L, "none"), list(60L, "none"))) {
    df <- .spg_df(pr, cfg[[1]], cfg[[2]])
    expect_length(df, ncol(pr$E))
  }
  # One vectorised call per fit, every voxel solved exactly (status 0), and
  # the inputs do clip eigenvalues.
  expect_equal(calls, 3L)
  expect_true(length(statuses) > 0L && all(statuses == 0L))
  expect_true(any(clipped > 0L))
})

test_that("the equivalence inputs reach every branch of the banded solver", {
  # path: 0 Cholesky certificate, 1 eigenvalues computed but nothing clipped,
  # 2 clipped-part correction, 3 rebuilt from the kept part (most clipped).
  # The values themselves are checked against the dense reference in
  # test-satterthwaite-psd-equivalence.R; this pins that those inputs cover
  # each branch, so a branch cannot silently go untested.
  paths <- integer(0)
  orig_cpp <- fmrireg:::satterthwaite_psd_traces_cpp
  testthat::local_mocked_bindings(
    satterthwaite_psd_traces_cpp = function(...) {
      res <- orig_cpp(...)
      paths <<- c(paths, res$path)
      res
    }
  )
  pr <- .spg_problem(n = 120L)
  pr$runs <- list(1:60, 61:120)
  for (cfg in list(list(30L, "none"), list(3L, "none"), list(8L, "tukey"))) {
    .spg_df(pr, cfg[[1]], cfg[[2]])
  }
  pr$E <- pr$E * 1e-4
  .spg_df(pr, 8L, "none")
  expect_true(all(c(0L, 2L, 3L) %in% paths))
})

test_that("a failed banded solve falls back to the dense per-block helper", {
  # The fallback triggers when LAPACK dstebz/dstein report info != 0, when
  # the Sturm count and the bisection disagree on the clipped count, or when
  # the autocovariances are not finite (C++ status != 0). Simulate a failure
  # for every voxel and check the result is unchanged.
  pr <- .spg_problem()
  expected <- .spg_df(pr)
  dense_calls <- 0L
  orig_cpp <- fmrireg:::satterthwaite_psd_traces_cpp
  orig_dense <- fmrireg:::.fmri_lm_projected_traces_dense
  testthat::local_mocked_bindings(
    satterthwaite_psd_traces_cpp = function(...) {
      res <- orig_cpp(...)
      res$status[] <- 1L
      res$tr_RV[] <- NA_real_
      res$tr_RVRV[] <- NA_real_
      res
    },
    .fmri_lm_projected_traces_dense = function(...) {
      dense_calls <<- dense_calls + 1L
      orig_dense(...)
    }
  )
  got <- .spg_df(pr)
  expect_equal(dense_calls, .spg_n_projected(pr))
  expect_lt(max(abs(got - expected) / expected), 1e-10)
})
