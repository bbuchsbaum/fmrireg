# Pass 3, stream R2-1: the voxelwise-covariate OLS path (ols_t_vcov_cpp)
# drops non-finite subjects per voxel, like ols_t_cpp, with per-voxel df.

make_block_gd_r21 <- function(Y, covars) {
  structure(
    list(blocks = list(list(Y = Y, V = NULL, covars = covars))),
    class = "group_data"
  )
}

r21_data <- function() {
  set.seed(211)
  S <- 15L
  P <- 5L
  covars <- data.frame(age = stats::rnorm(S, 40, 8))
  C <- matrix(stats::rnorm(S * P), S, P)
  Y <- matrix(stats::rnorm(S * P), S, P) + 0.5 * C + 0.02 * covars$age
  Y[3, 2] <- NA            # missing outcome
  C[c(1, 7), 3] <- c(NA, Inf)  # missing covariate
  Y[10, 4] <- NaN
  C[11, 4] <- NA           # both
  list(S = S, P = P, covars = covars, C = C, Y = Y)
}

test_that("R2-1: ols_t_vcov_cpp matches lm on complete cases per voxel", {
  d <- r21_data()
  X <- cbind(1, d$covars$age)
  ols <- ols_t_vcov_cpp(d$Y, X, d$C)
  expect_length(ols$df, d$P)
  for (j in seq_len(d$P)) {
    keep <- is.finite(d$Y[, j]) & is.finite(d$C[, j])
    m <- stats::lm(d$Y[keep, j] ~ d$covars$age[keep] + d$C[keep, j])
    cf <- summary(m)$coefficients
    expect_equal(unname(ols$beta[, j]), unname(cf[, "Estimate"]), tolerance = 1e-10, info = j)
    expect_equal(unname(ols$se[, j]), unname(cf[, "Std. Error"]), tolerance = 1e-10, info = j)
    expect_equal(unname(ols$t[, j]), unname(cf[, "t value"]), tolerance = 1e-10, info = j)
    expect_equal(ols$df[j], m$df.residual, info = j)
    expect_equal(ols$n_obs[j], sum(keep), info = j)
  }
})

test_that("R2-1: ols_t_vcov_cpp complete data keeps scalar df and lm results", {
  set.seed(212)
  S <- 12L
  Y <- matrix(stats::rnorm(S * 3), S, 3)
  C <- matrix(stats::rnorm(S * 3), S, 3)
  X <- cbind(1, stats::rnorm(S))
  ols <- ols_t_vcov_cpp(Y, X, C)
  expect_length(ols$df, 1L)
  expect_equal(ols$df, S - 3)
  for (j in 1:3) {
    m <- stats::lm(Y[, j] ~ X[, 2] + C[, j])
    expect_equal(unname(ols$beta[, j]), unname(stats::coef(m)), tolerance = 1e-10)
  }
})

test_that("R2-1: voxel left with too few or collinear subjects returns NaN", {
  set.seed(213)
  S <- 6L
  Y <- matrix(stats::rnorm(S * 2), S, 2)
  C <- matrix(stats::rnorm(S * 2), S, 2)
  X <- cbind(1, stats::rnorm(S))
  Y[1:3, 1] <- NA   # 3 left for 3 parameters -> no df
  ols <- ols_t_vcov_cpp(Y, X, C)
  expect_true(all(is.nan(ols$beta[, 1])))
  expect_true(is.nan(ols$df[1]))
  expect_true(all(is.finite(ols$beta[, 2])))
})

test_that("R2-1: fmri_ttest with voxelwise_cov containing NA matches lm per voxel", {
  d <- r21_data()
  gd <- make_block_gd_r21(d$Y, d$covars)
  fit <- fmri_ttest(gd, formula = ~ 1 + age, engine = "classic",
                    voxelwise_cov = d$C, voxel_name = "cov")
  expect_equal(rownames(fit$beta), c("(Intercept)", "age", "cov"))
  for (j in seq_len(d$P)) {
    keep <- is.finite(d$Y[, j]) & is.finite(d$C[, j])
    cc <- d$C[keep, j] - mean(d$C[keep, j])  # centred over the subjects used
    age <- d$covars$age[keep]
    m <- stats::lm(d$Y[keep, j] ~ age + cc)
    cf <- summary(m)$coefficients
    expect_equal(unname(fit$beta[, j]), unname(cf[, "Estimate"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$se[, j]), unname(cf[, "Std. Error"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$p[, j]), unname(cf[, "Pr(>|t|)"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df[, j]), rep(m$df.residual, 3), info = j)
  }
})

test_that("R2-1: paired test with NA in voxelwise covariate matches lm per voxel", {
  d <- r21_data()
  gd <- make_block_gd_r21(d$Y, data.frame(subject = seq_len(d$S)))
  fit <- suppressMessages(fmri_ttest(gd, formula = ~ 1, paired = TRUE, engine = "classic",
                                     voxelwise_cov = d$C, voxel_name = "cov"))
  for (j in seq_len(d$P)) {
    keep <- is.finite(d$Y[, j]) & is.finite(d$C[, j])
    cc <- d$C[keep, j] - mean(d$C[keep, j])
    m <- stats::lm(d$Y[keep, j] ~ cc)
    cf <- summary(m)$coefficients
    expect_equal(unname(fit$beta[, j]), unname(cf[, "Estimate"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$t[, j]), unname(cf[, "t value"]), tolerance = 1e-10, info = j)
    expect_equal(unname(fit$df[, j]), rep(m$df.residual, 2), info = j)
  }
})
