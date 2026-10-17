# Stream B of the 2026-10 simplification: inference helpers.
# Every test compares against a closed form, `lm()`/`anova()`, or the
# scale-invariance of the t statistic.

# B10 ---------------------------------------------------------------------

test_that("B10: voxelwise beta statistics are scale invariant and match lm()", {
  # beta_stats_matrix_voxelwise is the voxelwise-AR path's beta inference.
  # Its old absolute threshold (SE < sqrt(eps) -> t = 0) zeroed every t for
  # Y * 1e-8. Reference: summary(lm()) per voxel, at both scales.
  set.seed(20261017)
  n <- 50L
  X <- cbind(1, x1 = rnorm(n), x2 = rnorm(n))
  colnames(X) <- c("int", "x1", "x2")
  Y <- matrix(rnorm(n * 3L), n, 3L) + drop(X %*% c(1, 0.4, -0.3))
  stats_at <- function(Ys) {
    proj <- fmrireg:::.fast_preproject(X)
    B <- qr.coef(proj$qr, Ys)
    sigma <- sqrt(colSums(qr.resid(proj$qr, Ys)^2) / proj$dfres)
    fmrireg:::beta_stats_matrix_voxelwise(
      B, rep(list(proj$XtXinv), ncol(Ys)), sigma, proj$dfres, colnames(X)
    )$data[[1L]]
  }
  lm_t <- function(Ys) {
    t(vapply(seq_len(ncol(Ys)), function(v) {
      summary(stats::lm(Ys[, v] ~ X - 1))$coefficients[, "t value"]
    }, numeric(3)))
  }
  ref <- stats_at(Y)
  for (s in c(1e-8, 1e-12)) {
    small <- stats_at(Y * s)
    expect_equal(unname(small$stat[[1L]]), unname(lm_t(Y * s)), tolerance = 1e-8)
    expect_equal(small$stat[[1L]], ref$stat[[1L]], tolerance = 1e-8)
    expect_equal(small$prob[[1L]], ref$prob[[1L]], tolerance = 1e-8)
  }
})

test_that("B10: voxel status classification does not depend on data scale", {
  set.seed(4)
  y <- rnorm(30)
  for (s in c(1, 1e-8, 1e8)) {
    expect_identical(fmrireg:::.fmri_lm_voxel_status(y * s), "ok")
    expect_identical(fmrireg:::.fmri_lm_voxel_status(rep(3, 30) * s), "constant")
  }
  expect_identical(fmrireg:::.fmri_lm_voxel_status(numeric(30)), "all_zero")
  # A large offset with a resolvable (if small) fluctuation is still data.
  expect_identical(fmrireg:::.fmri_lm_voxel_status(1e8 + y * 0.1), "ok")
})

test_that("B10: an exactly degenerate standard error is handled alike on all paths", {
  # Coefficient 2 has zero (co)variance: its SE is exactly 0. `.safe_tstat`
  # defines t = 0 there (p = 1); every per-voxel-covariance path must agree.
  p <- 2L
  V <- 3L
  B <- rbind(c(1, -2, 0.5), c(0.3, 0.7, -1.1))   # p x V
  cv <- diag(c(0.25, 0))
  sigma <- rep(1, V)
  df <- 10
  expected_t <- rbind(B[1L, ] / 0.5, 0)           # p x V
  expected_p <- 2 * stats::pt(-abs(expected_t), df)

  # 1. beta_stats_matrix_voxelwise (per-voxel XtXinv)
  bs <- fmrireg:::beta_stats_matrix_voxelwise(
    B, rep(list(cv), V), sigma, df, c("a", "b")
  )$data[[1L]]
  expect_equal(unname(bs$stat[[1L]]), t(expected_t))
  expect_equal(unname(bs$prob[[1L]]), t(expected_p))

  # 2. fit_lm_contrasts_voxelwise
  l <- c(0, 1); attr(l, "colind") <- 1:2
  cres <- fmrireg:::fit_lm_contrasts_voxelwise(
    B, sigma^2, rep(list(cv), V), list(cb = l), list(), df
  )$cb$data[[1L]]
  expect_equal(cres$stat, rep(0, V))
  expect_equal(cres$prob, rep(1, V))

  # 3. .fmri_lm_update_inference (sandwich/HAC covariance)
  result <- list(
    betas = tibble::tibble(
      type = "beta", name = "parameter_estimates",
      data = list(tibble::tibble(
        estimate = list(t(B)), se = list(matrix(1, V, p)),
        stat = list(matrix(0, V, p)), prob = list(matrix(1, V, p)),
        sigma = list(sigma)
      ))
    ),
    contrasts = tibble::tibble(
      type = "contrast", name = "cb", conmat = list(c(0, 1)),
      colind = list(1:2),
      data = list(tibble::tibble(estimate = numeric(V), se = numeric(V),
                                 stat = numeric(V), prob = numeric(V)))
    )
  )
  upd <- fmrireg:::.fmri_lm_update_inference(result, rep(list(cv), V), df)
  ub <- upd$betas$data[[1L]]
  expect_equal(ub$stat[[1L]], t(expected_t))
  expect_equal(ub$prob[[1L]], t(expected_p))
  expect_equal(upd$contrasts$data[[1L]]$stat, rep(0, V))
  expect_equal(upd$contrasts$data[[1L]]$prob, rep(1, V))

  # 4. fit_contrasts.fmri_lm with a structured (per-voxel) covariance
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  bnames <- colnames(coef(fit))
  Vd <- nrow(as.matrix(coef(fit)))
  fit$result$variance_model <- structure(
    list(df_inference = rep(df, Vd), covariance_scope = "voxel",
         covariance = rep(list(cv), Vd)),
    class = "fmri_lm_variance_model"
  )
  w <- stats::setNames(c(0, 1), bnames)
  out <- fit_contrasts(fit, list(cb = w))$cb
  expect_equal(out$stat, rep(0, Vd))
  expect_equal(out$prob, rep(1, Vd))
  # ... and a non-degenerate contrast is the closed form est / se.
  w1 <- stats::setNames(c(1, 0), bnames)
  out1 <- fit_contrasts(fit, list(ca = w1))$ca
  est1 <- as.matrix(coef(fit))[, 1L]
  expect_equal(out1$stat, unname(est1 / 0.5))
  expect_equal(out1$prob, unname(2 * stats::pt(-abs(est1 / 0.5), df)))
})

# Contrast engine refactors (pinned against lm()/anova()) -------------------

.simplify_B_lm_fixture <- function(seed = 31L) {
  set.seed(seed)
  n <- 40L
  X <- cbind(1, a = rnorm(n), b = rnorm(n), c = rnorm(n))
  colnames(X) <- c("int", "a", "b", "c")
  Y <- X %*% matrix(rnorm(4L * 3L), 4L, 3L) + matrix(rnorm(n * 3L), n, 3L)
  list(X = X, Y = Y, n = n, p = 4L, V = 3L)
}

.simplify_B_lm_reference <- function(f) {
  t_ab <- f_bc <- numeric(f$V)
  p_t <- p_F <- numeric(f$V)
  for (v in seq_len(f$V)) {
    d <- data.frame(y = f$Y[, v], a = f$X[, "a"], b = f$X[, "b"], c = f$X[, "c"])
    # t for a - b: refit with the contrast as a single coefficient.
    full <- stats::lm(y ~ I(a) + I(a + b) + c, data = d)
    t_ab[v] <- summary(full)$coefficients["I(a)", "t value"]
    p_t[v] <- summary(full)$coefficients["I(a)", "Pr(>|t|)"]
    # F for b = c = 0 against the reduced model.
    av <- stats::anova(stats::lm(y ~ a, data = d), stats::lm(y ~ a + b + c, data = d))
    f_bc[v] <- av$F[2L]
    p_F[v] <- av$`Pr(>F)`[2L]
  }
  list(t = t_ab, p_t = p_t, F = f_bc, p_F = p_F)
}

test_that("compute_lm_contrasts matches lm()/anova() for t and F contrasts", {
  f <- .simplify_B_lm_fixture()
  ref <- .simplify_B_lm_reference(f)
  qx <- qr(f$X)
  B <- qr.coef(qx, f$Y)
  s2 <- colSums(qr.resid(qx, f$Y)^2) / (f$n - f$p)
  Fm <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1))
  out <- compute_lm_contrasts(
    B, solve(crossprod(f$X)), df = f$n - f$p, sigma2 = s2,
    t_contrasts = list(a_b = c(a = 1, b = -1)),
    f_contrasts = list(bc = Fm)
  )
  tt <- out$data[[which(out$name == "a_b")]]
  FF <- out$data[[which(out$name == "bc")]]
  expect_equal(tt$stat, ref$t, tolerance = 1e-8)
  expect_equal(tt$prob, ref$p_t, tolerance = 1e-8)
  expect_equal(FF$stat, ref$F, tolerance = 1e-8)
  expect_equal(FF$prob, ref$p_F, tolerance = 1e-8)
})

test_that("compute_lm_contrasts_from_suffstats equals compute_lm_contrasts and lm()", {
  f <- .simplify_B_lm_fixture(32L)
  ref <- .simplify_B_lm_reference(f)
  XtX <- crossprod(f$X)
  XtS <- crossprod(f$X, f$Y)
  StS <- colSums(f$Y^2)
  Fm <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1))
  args <- list(t_contrasts = list(a_b = c(a = 1, b = -1)),
               f_contrasts = list(bc = Fm))
  ss <- do.call(compute_lm_contrasts_from_suffstats,
                c(list(XtX = XtX, XtS = XtS, StS = StS, df = f$n - f$p), args))
  tt <- ss$data[[which(ss$name == "a_b")]]
  FF <- ss$data[[which(ss$name == "bc")]]
  expect_equal(tt$stat, ref$t, tolerance = 1e-8)
  expect_equal(tt$prob, ref$p_t, tolerance = 1e-8)
  expect_equal(FF$stat, ref$F, tolerance = 1e-8)
  expect_equal(FF$prob, ref$p_F, tolerance = 1e-8)

  XtXinv <- chol2inv(chol(XtX))
  B <- XtXinv %*% XtS
  s2 <- pmax(pmax(StS - colSums(B * XtS), 0) / (f$n - f$p), .Machine$double.eps)
  direct <- do.call(compute_lm_contrasts,
                    c(list(B = B, XtXinv = XtXinv, df = f$n - f$p, sigma2 = s2,
                           columns = colnames(f$X)), args))
  expect_equal(ss, direct)

  # Supplied sigma takes precedence over the suff-stat residual variance.
  ss_sigma <- do.call(compute_lm_contrasts_from_suffstats,
                      c(list(XtX = XtX, XtS = XtS, StS = StS, df = f$n - f$p,
                             sigma = rep(2, f$V)), args))
  direct_sigma <- do.call(compute_lm_contrasts,
                          c(list(B = B, XtXinv = XtXinv, df = f$n - f$p,
                                 sigma = rep(2, f$V), columns = colnames(f$X)), args))
  expect_equal(ss_sigma, direct_sigma)
})

test_that("reference lm-object contrast helpers match summary.lm()", {
  f <- .simplify_B_lm_fixture(33L)
  fit <- stats::lm(f$Y[, 1] ~ f$X - 1)
  cu <- fmrireg:::.lm_cov_unscaled(fit, varnames = colnames(f$X))
  expect_equal(matrix(cu, 4L), unname(summary(fit)$cov.unscaled), tolerance = 1e-10)
  bs <- fmrireg:::beta_stats(fit, colnames(f$X))
  expect_equal(unname(drop(bs$data[[1L]]$stat[[1L]])),
               unname(summary(fit)$coefficients[, "t value"]), tolerance = 1e-10)
})

# Exact result-field matching ----------------------------------------------

test_that("fit_contrasts.fmri_lm does not read sigma2 as sigma", {
  # `result$sigma` partial-matches `sigma2` when `sigma` is absent, so the
  # old code silently used the variance as the standard deviation. With
  # exact matching a result without `sigma` (the beta table carries no
  # sigma column) fails loudly instead.
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  bnames <- colnames(coef(fit))
  w <- stats::setNames(c(1, -1), bnames)
  alt <- fit
  sigma <- as.numeric(alt$result[["sigma"]])
  alt$result[["sigma"]] <- NULL
  alt$result[["sigma2"]] <- (3 * sigma)^2
  expect_error(fit_contrasts(alt, list(d = w)), "Length of sigma")

  # With `sigma` present, the stored scale is used: se = sigma * sqrt(l'XtXinv l).
  ref <- fit_contrasts(fit, list(d = w))$d
  XtXinv <- as.matrix(fit$result[["cov.unscaled"]])
  ei <- fit$result[["event_indices"]]
  if (nrow(XtXinv) != length(ei)) XtXinv <- XtXinv[ei, ei, drop = FALSE]
  l <- c(1, -1)
  expect_equal(ref$se, sigma * sqrt(drop(crossprod(l, XtXinv %*% l))))
})

test_that("tidy.fmri_lm does not read df.residual as the df record", {
  # `result$df` partial-matched `df.residual` when `df` was absent and then
  # failed on `$inference` of a number.
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  ref <- tidy(fit)
  alt <- fit
  alt$result[["df"]] <- NULL
  alt$result[["df.residual"]] <- 5
  got <- tidy(alt)
  expect_equal(got$estimate, ref$estimate)
  expect_equal(got$statistic, ref$statistic)
  expect_false("df_inference" %in% names(got))
})
