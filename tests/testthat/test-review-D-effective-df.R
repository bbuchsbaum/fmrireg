# Review 2026-10, stream D, item D4: one effective-df implementation.
#
# Robust inference uses df = sum(w) - rank (CLAUDE.md, "Effective df for
# robust models"): a weight-w observation contributes w of an observation, so
# a zero weight contributes nothing and 0/1 weights reproduce lm.wfit()'s
# residual df exactly.

test_that("robust effective df is sum(w) - rank (hand computation)", {
  w <- c(rep(1, 90), rep(0, 10))
  # 0/1 weights: exactly the WLS residual df that lm() reports.
  set.seed(1)
  X <- cbind(1, rnorm(100), rnorm(100), rnorm(100), rnorm(100))
  y <- rnorm(100)
  ref <- lm(y ~ X - 1, weights = w)
  expect_equal(fmrireg:::calculate_effective_df(100, 5, robust_weights = w),
               ref$df.residual)
  expect_equal(fmrireg:::calculate_effective_df(100, 5, robust_weights = w), 85)

  w2 <- c(rep(1, 60), rep(0.5, 30), rep(0.2, 10))
  expect_equal(fmrireg:::calculate_effective_df(100, 5, robust_weights = w2),
               sum(w2) - 5)
  # n is implied by the weights; it is not double-counted.
  expect_equal(fmrireg:::calculate_effective_df(73, 5, robust_weights = w2),
               sum(w2) - 5)
  # All-ones weights are a no-op.
  expect_equal(fmrireg:::calculate_effective_df(100, 5, robust_weights = rep(1, 100)), 95)
  expect_equal(fmrireg:::.robust_effective_df(w2, rank = 5), sum(w2) - 5)
})

test_that("robust and Satterthwaite adjustments compose by the robust ratio", {
  set.seed(15)
  n <- 120L
  X <- cbind(1, rnorm(n), rnorm(n))
  V <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  w <- c(rep(1, 100), rep(0.5, 20))
  base <- fmrireg:::calculate_effective_df(n, 3, method = "satterthwaite",
                                           X = X, resid_cov = V)
  both <- fmrireg:::calculate_effective_df(n, 3, robust_weights = w,
                                           method = "satterthwaite",
                                           X = X, resid_cov = V)
  expect_equal(both, base * (sum(w) - 3) / (n - 3), tolerance = 1e-10)
})

test_that("robust fits report sum(w) - rank as their inference df", {
  set.seed(7)
  n <- 100L
  ev <- data.frame(onset = seq(6, 170, by = 14),
                   condition = factor(rep(c("A", "B"), length.out = 12)),
                   run = 1)
  Y <- matrix(rnorm(n * 4, mean = 10), n, 4)
  Y[c(5, 30, 31), ] <- Y[c(5, 30, 31), ] + 30
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  for (scope in c("joint", "runwise_meta")) {
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                   durations = 0,
                   control = fmri_lm_control(
                     estimation = estimation_spec(scope),
                     robust = robust_spec("bisquare")
                   ),
                   compute = compute_spec(voxel_chunks = 2))
    w <- as.numeric(unlist(fit$result$robust_weights))
    rank <- qr(as.matrix(design_matrix(fit$model)))$rank
    expect_equal(unique(fit$result$df$inference), sum(w) - rank, info = scope)
    # Beta p-values use the same df.
    bt <- fit$result$betas$data[[1]]
    expect_equal(c(bt$prob[[1]]),
                 c(2 * pt(-abs(bt$stat[[1]]), sum(w) - rank)),
                 tolerance = 1e-10, info = scope)
  }
})

test_that("unreachable effective-df copies are gone", {
  ns <- asNamespace("fmrireg")
  for (fn in c("compute_effective_df", "compute_effective_df_ar",
               "compute_effective_df_robust", "compute_effective_df_ar_robust",
               "compute_robust_effective_df", "residual_effective_df",
               "satterthwaite_df", "kenward_roger_df",
               "adjust_stats_for_effective_df")) {
    expect_false(exists(fn, envir = ns, inherits = FALSE), info = fn)
  }
})
