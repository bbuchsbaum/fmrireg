library(fmrireg)
library(testthat)

test_that("full-rank near-collinear design matches QR reference", {
  n <- 120
  set.seed(2026)

  x1 <- as.numeric(scale(rnorm(n)))
  z <- rnorm(n)
  z <- z - x1 * sum(x1 * z) / sum(x1^2)
  x2 <- x1 + 3e-7 * z
  X <- cbind(1, x1, x2)

  expect_equal(qr(X)$rank, ncol(X))
  expect_gt(kappa(X), 1e6)

  Y <- matrix(rnorm(n), ncol = 1)

  proj <- fmrireg:::.fast_preproject(X)
  ctx <- fmrireg:::glm_context(X = X, Y = Y, proj = proj)
  fit <- fmrireg:::solve_glm_core(ctx)

  beta_ref <- qr.solve(X, Y)
  rel_err <- max(abs(fit$betas - beta_ref)) / max(1, max(abs(beta_ref)))

  expect_lt(rel_err, 1e-6)
})
