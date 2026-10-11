# Review 2026-10, stream E7: residual df after soft subspace projection.
# soft_projection() computes effective_df but the fitting code keeps
# dfres = n - p. .soft_projection_residual_df() gives the unbiased value; the
# wiring into the runwise fitter belongs to fmri_lm_strategies.R (stream D).

test_that("E7: residual df equals tr((I - H) P^2) computed densely", {
  set.seed(1)
  n <- 80; k <- 12
  N <- matrix(rnorm(n * k), n, k)
  X <- cbind(1, rnorm(n), sin(seq_len(n) / 5))
  for (lam in list("auto", 0.5, 1e3)) {
    proj <- soft_projection(N, lambda = lam)
    Xp <- proj$P_lambda(X)
    P <- proj$P_lambda(diag(n))
    H <- Xp %*% solve(crossprod(Xp), t(Xp))
    dense <- sum(diag((diag(n) - H) %*% P %*% P))
    expect_equal(fmrireg:::.soft_projection_residual_df(proj, Xp), dense,
                 tolerance = 1e-8)
    # The naive n - p overstates the df whenever the projection removes anything.
    expect_lt(dense, n - ncol(X))
  }
})

test_that("E7: residual df gives an unbiased sigma^2 (Monte Carlo)", {
  set.seed(2)
  n <- 60; k <- 10; nrep <- 4000
  N <- matrix(rnorm(n * k), n, k)
  X <- cbind(1, rnorm(n))
  proj <- soft_projection(N, lambda = "auto")
  Xp <- proj$P_lambda(X)
  E <- matrix(rnorm(n * nrep), n, nrep)          # sigma = 1
  Ep <- proj$P_lambda(E)
  rss <- colSums(qr.resid(qr(Xp), Ep)^2)
  df_new <- fmrireg:::.soft_projection_residual_df(proj, Xp)
  # mean(rss / df) ~ 1 with Monte Carlo SE ~ sqrt(2/df/nrep) ~ 0.004
  expect_equal(mean(rss / df_new), 1, tolerance = 0.02)
  # The df currently used by the fitter (n - p) is biased low by > 10%.
  expect_lt(mean(rss / (n - ncol(X))), 0.9)
})

test_that("E7: limits -- lambda -> Inf gives n - p, hard projection gives n - k - p", {
  set.seed(3)
  n <- 50; k <- 5
  N <- matrix(rnorm(n * k), n, k)
  X <- cbind(1, rnorm(n))
  p_inf <- soft_projection(N, lambda = 1e12)
  expect_equal(fmrireg:::.soft_projection_residual_df(p_inf, p_inf$P_lambda(X)),
               n - 2, tolerance = 1e-6)
  p_hard <- soft_projection(N, lambda = 1e-12)
  expect_equal(fmrireg:::.soft_projection_residual_df(p_hard, p_hard$P_lambda(X)),
               n - k - 2, tolerance = 1e-6)
})
