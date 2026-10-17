# Review 2026-10, stream C: contrast statistics.
#
# C7: absolute thresholds (se < sqrt(eps), |sigma2| < sqrt(eps)) zeroed valid
#     statistics for small-scale data.
# C8: `1 - pf()` instead of the upper tail.
# D5 (contrast part): effective df built from n <- dfres + ncol(X), which is
#     wrong when the design is rank deficient (dfres = n - rank).
# P-C2: reference F path looped over voxels and recomputed cov.unscaled.
# P-C3: compute_contrast() filled standard errors with a voxel loop.

.cs_design <- function(n = 60L, seed = 5L) {
  set.seed(seed)
  X <- cbind(1, rnorm(n), rnorm(n), rnorm(n))
  colnames(X) <- c("int", "a", "b", "c")
  B <- c(0.5, 1, -0.5, 0.25)
  Y <- drop(X %*% B) + matrix(rnorm(n * 5), n, 5)
  list(X = X, Y = Y)
}

.cs_fast <- function(X, Y) {
  qrx <- qr(X)
  B <- qr.coef(qrx, Y)
  res <- qr.resid(qrx, Y)
  df <- nrow(X) - qrx$rank
  sigma2 <- colSums(res^2) / df
  XtXinv <- chol2inv(qr.R(qrx))
  list(B = B, sigma2 = sigma2, XtXinv = XtXinv, df = df)
}

test_that("C7: fast t, F and beta statistics are invariant to the scale of Y", {
  d <- .cs_design()
  l <- matrix(c(0, 1, -1, 0), 1)
  L <- rbind(c(0, 1, 0, 0), c(0, 0, 1, 0))
  base <- .cs_fast(d$X, d$Y)
  t0 <- fmrireg:::.fast_t_contrast(base$B, base$sigma2, base$XtXinv, l, base$df)
  F0 <- fmrireg:::.fast_F_contrast(base$B, base$sigma2, base$XtXinv, L, base$df)
  b0 <- fmrireg:::beta_stats_matrix(base$B, base$XtXinv, sqrt(base$sigma2),
                                    base$df, colnames(d$X))
  for (s in c(1e-8, 1e-4, 1e4)) {
    f <- .cs_fast(d$X, d$Y * s)
    ts <- fmrireg:::.fast_t_contrast(f$B, f$sigma2, f$XtXinv, l, f$df)
    Fs <- fmrireg:::.fast_F_contrast(f$B, f$sigma2, f$XtXinv, L, f$df)
    bs <- fmrireg:::beta_stats_matrix(f$B, f$XtXinv, sqrt(f$sigma2), f$df,
                                      colnames(d$X))
    expect_equal(ts$stat, t0$stat, tolerance = 1e-8, info = paste("t, scale", s))
    expect_equal(ts$prob, t0$prob, tolerance = 1e-8)
    expect_equal(Fs$stat, F0$stat, tolerance = 1e-8, info = paste("F, scale", s))
    expect_equal(Fs$prob, F0$prob, tolerance = 1e-8)
    expect_equal(bs$data[[1]]$stat[[1]], b0$data[[1]]$stat[[1]], tolerance = 1e-8,
                 info = paste("beta t, scale", s))
  }

  # Reference values from lm.
  m <- lm(d$Y[, 2] ~ 0 + d$X)
  tv <- summary(m)$coefficients[, "t value"]
  expect_equal(unname(b0$data[[1]]$stat[[1]][2, ]), unname(tv), tolerance = 1e-8)
})

test_that("C7/C8: reference lm path (beta_stats, fit_Fcontrasts) is scale invariant and matches anova", {
  d <- .cs_design()
  for (s in c(1, 1e-8, 1e-4, 1e4)) {
    Ys <- d$Y * s
    fit <- lm(Ys ~ 0 + d$X)
    bs <- fmrireg:::beta_stats(fit, colnames(d$X))
    Fc <- fmrireg:::fit_Fcontrasts(fit, diag(2), colind = 2:3)
    for (v in seq_len(ncol(Ys))) {
      mv <- lm(Ys[, v] ~ 0 + d$X)
      tv <- summary(mv)$coefficients[, "t value"]
      expect_equal(unname(bs$data[[1]]$stat[[1]][v, ]), unname(tv), tolerance = 1e-7,
                   info = paste("scale", s))
      a <- anova(lm(Ys[, v] ~ 0 + d$X[, c(1, 4)]), mv)
      expect_equal(Fc$stat[v], a$F[2], tolerance = 1e-7, info = paste("scale", s))
      expect_equal(Fc$prob[v], a$`Pr(>F)`[2], tolerance = 1e-7)
    }
  }
})

test_that("C8: F p-values keep precision in the far upper tail", {
  set.seed(9)
  n <- 200
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- drop(X %*% c(0, 0.6, 0.6)) + rnorm(n)
  fit <- lm(cbind(y, y) ~ 0 + X)
  Fc <- fmrireg:::fit_Fcontrasts(fit, diag(2), colind = 2:3)
  ref <- pf(Fc$stat, 2, n - 3, lower.tail = FALSE)
  # In this regime 1 - pf() rounds to exactly 0.
  expect_true(all(ref > 0 & ref < 1e-17))
  expect_equal(log(Fc$prob), log(ref))
})

test_that("C7: fmri_lm t and F statistics are invariant to rescaling the data", {
  set.seed(21)
  n <- 90L
  etab <- data.frame(onset = seq(4, 78, by = 6),
                     cond = factor(rep(c("A", "B", "C"), length.out = 13)),
                     run = 1L)
  Y <- matrix(rnorm(n * 3), n, 3)
  cons <- contrast_set(
    pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"),
    oneway_contrast(~ cond, name = "main")
  )
  fit_at <- function(s) {
    dset <- matrix_frame(Y * s, TR = 1, run_length = n, event_table = etab)
    fmri_lm(onset ~ hrf(cond, contrasts = cons), block = ~ run, dataset = dset)
  }
  f1 <- fit_at(1)
  for (s in c(1e-4, 1e4)) {
    fs <- fit_at(s)
    expect_equal(stats(fs, "contrasts"), stats(f1, "contrasts"), tolerance = 1e-6,
                 info = paste("scale", s))
    expect_equal(stats(fs, "F"), stats(f1, "F"), tolerance = 1e-6,
                 info = paste("scale", s))
    expect_equal(stats(fs, "betas"), stats(f1, "betas"), tolerance = 1e-6,
                 info = paste("scale", s))
  }
})

test_that("D5: effective-df inputs use the design rank, and rank-deficient df match lm", {
  d <- .cs_design()
  Xa <- cbind(d$X, alias = d$X[, "a"] + d$X[, "b"])   # aliased column
  qrx <- qr(Xa)
  expect_equal(qrx$rank, 4L)
  B <- matrix(0, ncol(Xa), ncol(d$Y))
  B[qrx$pivot[1:4], ] <- qr.coef(qrx, d$Y)[qrx$pivot[1:4], ]
  df <- nrow(Xa) - qrx$rank
  sigma2 <- colSums(qr.resid(qrx, d$Y)^2) / df
  XtXinv <- matrix(0, 5, 5)
  XtXinv[1:4, 1:4] <- chol2inv(qr.R(qrx)[1:4, 1:4])
  attr(XtXinv, "aliased") <- 5L
  attr(B, "aliased") <- 5L

  # The robust df subtracts the design rank (4), not its width (5).
  seen <- list()
  local_mocked_bindings(
    .robust_effective_df = function(weights, rank, df_base = NULL) {
      seen[[length(seen) + 1L]] <<- c(n = length(weights), p = rank)
      sum(weights) - rank
    }
  )
  w <- rep(1, nrow(Xa))
  l <- matrix(c(0, 1, -1, 0, 0), 1)
  tt <- fmrireg:::.fast_t_contrast(B, sigma2, XtXinv, l, df, robust_weights = w,
                                   ar_order = 1, aliased = 5L)
  FF <- fmrireg:::.fast_F_contrast(B, sigma2, XtXinv,
                                   rbind(c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0)), df,
                                   robust_weights = w, ar_order = 1, aliased = 5L)
  bs_w <- fmrireg:::beta_stats_matrix(B, XtXinv, sqrt(sigma2), df, colnames(Xa),
                                      robust_weights = w, ar_order = 1)
  expect_length(seen, 3L)
  for (s in seen) {
    expect_equal(unname(s["n"]), nrow(Xa))
    expect_equal(unname(s["p"]), qrx$rank)
  }
  bs <- fmrireg:::beta_stats_matrix(B, XtXinv, sqrt(sigma2), df, colnames(Xa),
                                    ar_order = 1)
  # Unit weights give exactly the residual df.
  expect_equal(bs_w$data[[1]]$prob[[1]], bs$data[[1]]$prob[[1]])

  # And the statistics agree with lm on the rank-deficient design.
  for (v in seq_len(ncol(d$Y))) {
    m <- lm(d$Y[, v] ~ 0 + Xa)
    expect_equal(m$df.residual, df)
    cf <- summary(m)$coefficients
    expect_equal(unname(bs$data[[1]]$stat[[1]][v, 1:4]), unname(cf[, "t value"]),
                 tolerance = 1e-8)
    expect_equal(unname(bs$data[[1]]$prob[[1]][v, 1:4]), unname(cf[, "Pr(>|t|)"]),
                 tolerance = 1e-8)
    est <- coef(m)[2] - coef(m)[3]
    se <- sqrt(drop(c(1, -1) %*% vcov(m)[2:3, 2:3] %*% c(1, -1)))
    expect_equal(tt$stat[v], unname(est / se), tolerance = 1e-8)
    expect_equal(tt$prob[v], unname(2 * pt(-abs(est / se), df)), tolerance = 1e-8)
  }
})

test_that("P-C2: reference F path equals .fast_F_contrast and anova(lm)", {
  d <- .cs_design(n = 80L, seed = 77L)
  fit <- lm(d$Y ~ 0 + d$X)
  # cell-oriented (k x r) weights, as from oneway_contrast()
  W <- rbind(c(1, 0), c(-1, 1), c(0, -1))
  Fc <- fmrireg:::fit_Fcontrasts(fit, W, colind = 2:4)
  f <- .cs_fast(d$X, d$Y)
  L <- matrix(0, 2, 4)
  L[, 2:4] <- t(W)
  Ff <- fmrireg:::.fast_F_contrast(f$B, f$sigma2, f$XtXinv, L, f$df)
  expect_equal(Fc$stat, Ff$stat, tolerance = 1e-10)
  expect_equal(Fc$prob, Ff$prob, tolerance = 1e-10)
  expect_equal(Fc$estimate, Ff$estimate, tolerance = 1e-10)
  for (v in seq_len(ncol(d$Y))) {
    full <- lm(d$Y[, v] ~ 0 + d$X)
    red <- lm(d$Y[, v] ~ 0 + d$X[, 1] + I(d$X[, 2] + d$X[, 3] + d$X[, 4]))
    a <- anova(red, full)
    expect_equal(Fc$stat[v], a$F[2], tolerance = 1e-8)
    expect_equal(Fc$prob[v], a$`Pr(>F)`[2], tolerance = 1e-8)
  }
  # A precomputed cov.unscaled gives identical results.
  cu <- fmrireg:::.lm_cov_unscaled(fit)
  Fc2 <- fmrireg:::fit_Fcontrasts(fit, W, colind = 2:4, cov.unscaled = cu)
  expect_identical(Fc2$stat, Fc$stat)
})

