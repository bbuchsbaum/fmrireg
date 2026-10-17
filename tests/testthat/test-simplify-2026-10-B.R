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
