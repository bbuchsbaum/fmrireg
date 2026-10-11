# Review 2026-10, stream B: classic t-test engine (B7) and paired voxelwise
# covariates (B8).

make_block_gd <- function(Y, covars, V = NULL) {
  structure(
    list(blocks = list(list(Y = Y, V = V, covars = covars))),
    class = "group_data"
  )
}

test_that("B7: classic one-sample t-test drops non-finite subjects per voxel (vs t.test)", {
  set.seed(7)
  S <- 12L
  P <- 5L
  Y <- matrix(stats::rnorm(S * P, mean = 0.4), S, P)
  Y[3, 2] <- NA
  Y[c(1, 8), 4] <- c(NaN, Inf)
  gd <- make_block_gd(Y, data.frame(subject = seq_len(S)))
  fit <- fmri_ttest(gd, formula = ~ 1, engine = "classic")

  for (j in seq_len(P)) {
    y <- Y[is.finite(Y[, j]), j]
    ref <- stats::t.test(y)
    expect_equal(unname(fit$beta[1, j]), mean(y), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$se[1, j]), stats::sd(y) / sqrt(length(y)), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$t[1, j]), unname(ref$statistic), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df[1, j]), unname(ref$parameter), info = j)
    expect_equal(unname(fit$p[1, j]), ref$p.value, tolerance = 1e-10, info = j)
  }
})

test_that("B7: classic two-sample/ANCOVA engine matches lm on complete cases, incl. contrasts", {
  set.seed(8)
  S <- 16L
  P <- 4L
  covars <- data.frame(
    group = factor(rep(c("A", "B"), each = S / 2)),
    age = stats::rnorm(S, 40, 8)
  )
  Y <- matrix(stats::rnorm(S * P), S, P) + 0.6 * (covars$group == "B")
  Y[2, 1] <- NA
  Y[c(5, 12, 13), 3] <- NA
  gd <- make_block_gd(Y, covars)
  fit <- fmri_ttest(gd, formula = ~ 1 + group + age, engine = "classic",
                    contrast = c(group = 1, age = 2))

  for (j in seq_len(P)) {
    keep <- is.finite(Y[, j])
    dat <- cbind(covars[keep, ], y = Y[keep, j])
    m <- stats::lm(y ~ 1 + group + age, data = dat)
    cf <- summary(m)$coefficients
    # default sign is A - B, raw lm coefficient is B - A
    expect_equal(unname(fit$beta["group", j]), -cf["groupB", "Estimate"], tolerance = 1e-10, info = j)
    expect_equal(unname(fit$se["group", j]), cf["groupB", "Std. Error"], tolerance = 1e-10, info = j)
    expect_equal(unname(fit$t["age", j]), cf["age", "t value"], tolerance = 1e-10, info = j)
    expect_equal(unname(fit$p["age", j]), cf["age", "Pr(>|t|)"], tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df["age", j]), m$df.residual, info = j)

    # contrast A - B + 2 * age, i.e. -groupB + 2 * age in raw coefficients
    L <- c(0, -1, 2)
    est <- sum(L * stats::coef(m))
    se <- sqrt(drop(t(L) %*% stats::vcov(m) %*% L))
    expect_equal(fit$beta_contrast[j], est, tolerance = 1e-10, info = j)
    expect_equal(fit$se_contrast[j], se, tolerance = 1e-10, info = j)
    expect_equal(fit$p_contrast[j],
                 2 * stats::pt(abs(est / se), m$df.residual, lower.tail = FALSE),
                 tolerance = 1e-10, info = j)
  }
})

test_that("B7: ols_t_cpp complete-data output is unchanged (scalar df)", {
  set.seed(9)
  Y <- matrix(stats::rnorm(30), 10, 3)
  X <- cbind(1, stats::rnorm(10))
  ols <- ols_t_cpp(Y, X)
  expect_length(ols$df, 1L)
  expect_equal(ols$df, 8)
  expect_equal(ols$beta, solve(crossprod(X), crossprod(X, Y)), tolerance = 1e-10)
})

test_that("B8: paired test with voxelwise covariate runs and matches lm per voxel", {
  set.seed(10)
  S <- 14L
  P <- 3L
  D <- matrix(stats::rnorm(S * P, mean = 0.3), S, P)
  C <- matrix(stats::rnorm(S * P), S, P)
  D <- D + 0.5 * C
  gd <- make_block_gd(D, data.frame(subject = seq_len(S)))
  fit <- suppressMessages(fmri_ttest(gd, formula = ~ 1, paired = TRUE, engine = "classic",
                                     voxelwise_cov = C, voxel_name = "cov"))
  expect_equal(rownames(fit$beta), c("(Intercept)", "cov"))
  for (j in seq_len(P)) {
    cc <- C[, j] - mean(C[, j])
    m <- stats::lm(D[, j] ~ cc)
    cf <- summary(m)$coefficients
    expect_equal(unname(fit$beta[, j]), unname(cf[, "Estimate"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$se[, j]), unname(cf[, "Std. Error"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$t[, j]), unname(cf[, "t value"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$p[, j]), unname(cf[, "Pr(>|t|)"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df[, j]), rep(m$df.residual, 2))
  }
})
