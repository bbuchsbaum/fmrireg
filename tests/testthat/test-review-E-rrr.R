# Review 2026-10, stream E: rrr_gls engine inference and RNG hygiene.

.review_e_rrr_data <- function(V = 2000L, seed = 1L) {
  set.seed(seed)
  ons <- seq(10, 380, 15)
  ev <- data.frame(
    onset = ons,
    condition = factor(rep(c("a", "b", "c"), length.out = length(ons))),
    run = 1
  )
  Y <- matrix(rnorm(200 * V), 200, V)
  list(
    ev = ev,
    Y = Y,
    dset = matrix_frame(Y, TR = 2, run_length = 200, event_table = ev)
  )
}

test_that("E1: full-rank rrr_gls betas and SEs equal OLS fmri_lm", {
  d <- .review_e_rrr_data(V = 300L)
  f_ols <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$dset)
  f_rrr <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$dset,
                   engine = "rrr_gls", engine_args = list(bootstrap_n = 0L))

  expect_equal(as.matrix(coef(f_rrr)), as.matrix(coef(f_ols)), tolerance = 1e-8)
  expect_equal(as.matrix(standard_error(f_rrr)), as.matrix(standard_error(f_ols)),
               tolerance = 1e-8)
  expect_equal(as.matrix(stats(f_rrr)), as.matrix(stats(f_ols)), tolerance = 1e-8)

  # Independent reference: lm() on voxel 1 with the full design.
  X <- as.matrix(design_matrix(f_ols$model))
  lmfit <- summary(stats::lm(d$Y[, 7] ~ X - 1))
  ev_cols <- seq_len(3)
  expect_equal(unname(as.matrix(standard_error(f_rrr))[7, ]),
               unname(lmfit$coefficients[ev_cols, "Std. Error"]), tolerance = 1e-8)
})

test_that("E1: reduced-rank rrr_gls t-statistics are calibrated on pure noise", {
  d <- .review_e_rrr_data(V = 2000L, seed = 11L)
  for (r in 1:2) {
    cset <- contrast_set(pair_contrast(~ condition == "a", ~ condition == "b", name = "ab"))
    fit <- fmri_lm(onset ~ hrf(condition, contrasts = cset), block = ~run,
                   dataset = d$dset, engine = "rrr_gls",
                   engine_args = list(rank = r))
    tt <- as.matrix(stats(fit))
    # Nominal 5%; a data-chosen reduced-rank subspace shrinks the estimate so
    # the test may be conservative, but must never be anti-conservative.
    # (Before the fix: ~90% of null voxels had |t| > 1.96, median |t| ~ 27.)
    expect_lt(mean(abs(tt) > 1.96), 0.08)
    expect_lt(stats::median(abs(tt)), 1)
    ct <- as.matrix(stats(fit, "contrasts"))
    expect_equal(ncol(ct), 1L)
    expect_lt(mean(abs(ct) > 1.96), 0.08)
  }

  # Full rank on the same null data is the OLS test: rate inside a binomial band.
  fit_full <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$dset,
                      engine = "rrr_gls")
  rate_full <- mean(abs(as.matrix(stats(fit_full))) > 1.96)
  expect_gt(rate_full, 0.035)
  expect_lt(rate_full, 0.065)
})

test_that("E3: rrr_gls bootstrap_seed does not alter the global RNG state", {
  d <- .review_e_rrr_data(V = 50L)
  set.seed(99)
  before <- .Random.seed
  fit1 <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$dset,
                  engine = "rrr_gls",
                  engine_args = list(rank = 2L, se_mode = "bootstrap",
                                     bootstrap_n = 10L, bootstrap_seed = 42L))
  expect_true(identical(.Random.seed, before))

  # Seeded bootstrap stays reproducible.
  fit2 <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$dset,
                  engine = "rrr_gls",
                  engine_args = list(rank = 2L, se_mode = "bootstrap",
                                     bootstrap_n = 10L, bootstrap_seed = 42L))
  expect_identical(as.matrix(standard_error(fit1)), as.matrix(standard_error(fit2)))
})

test_that("P-E1: bootstrap with hoisted design QR matches the naive refit loop", {
  set.seed(5)
  n <- 80L; k <- 3L; V <- 40L
  X0 <- matrix(rnorm(n * k), n, k)
  Y0 <- X0 %*% matrix(rnorm(k * V), k, V) + matrix(rnorm(n * V), n, V)
  rc <- list(rank_mode = "fixed", rank = 2L, energy_keep = 0.99, rss_budget = NULL)
  fit0 <- fmrireg:::.rrr_fit_task_subspace(X0, Y0, rc)
  cons <- list(ab = structure(c(1, -1), colind = 1:2))

  out <- fmrireg:::.rrr_bootstrap_task_se(
    X0, Y0, fit0$B_task, rc, cons, nboot = 25L, block_size = 4L, seed = 3L
  )

  # Naive reference: refit the full task subspace on each resampled dataset.
  naive <- local({
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
    set.seed(3L)
    fitted <- X0 %*% fit0$B_task
    resid <- Y0 - fitted
    starts <- seq.int(1L, n, by = 4L)
    blocks <- lapply(starts, function(s) seq.int(s, min(n, s + 3L)))
    Bs <- vector("list", 25L)
    for (b in 1:25) {
      idx <- integer(0)
      while (length(idx) < n) idx <- c(idx, blocks[[sample.int(length(blocks), 1L)]])
      idx <- idx[seq_len(n)]
      Bs[[b]] <- fmrireg:::.rrr_fit_task_subspace(X0, fitted + resid[idx, , drop = FALSE], rc)$B_task
    }
    arr <- simplify2array(Bs)
    se <- apply(arr, c(1, 2), stats::sd)
    tcon <- sapply(Bs, function(B) as.numeric(c(1, -1, 0) %*% B))
    list(se_beta = se, se_t = apply(tcon, 1, stats::sd))
  })

  expect_equal(out$se_beta, naive$se_beta, tolerance = 1e-10)
  expect_equal(out$se_t$ab, naive$se_t, tolerance = 1e-10)
})
