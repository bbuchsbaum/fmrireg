# Tests for fmri_ttest functionality

test_that("fmri_ttest works for one-sample t-test", {
  set.seed(123)
  
  # Simulate data
  S <- 20  # subjects
  P <- 100  # features
  Y <- matrix(rnorm(S * P, mean = 0.5), nrow = S, ncol = P)
  
  # Create simple group_data-like structure
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = NULL,
        covars = data.frame(subject = 1:S),
        feature = list(group = rep(1:10, each = 10))
      )
    )
  )
  class(gd) <- "group_data"
  
  # One-sample t-test
  fit <- fmri_ttest(gd, formula = ~ 1, engine = "classic")
  
  expect_s3_class(fit, "fmri_ttest_fit")
  expect_equal(dim(fit$beta), c(1, P))
  expect_equal(dim(fit$t), c(1, P))
  expect_equal(dim(fit$p), c(1, P))
  expect_equal(fit$n_subjects, S)
  expect_equal(fit$n_features, P)
})

test_that("fmri_ttest works for two-sample t-test", {
  set.seed(124)
  
  # Simulate data with group difference
  S <- 30
  P <- 50
  group <- factor(rep(c("A", "B"), each = 15))
  Y <- matrix(rnorm(S * P), nrow = S, ncol = P)
  Y[group == "B", ] <- Y[group == "B", ] + 0.5  # Add group effect
  
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = NULL,
        covars = data.frame(subject = 1:S, group = group),
        feature = NULL
      )
    )
  )
  class(gd) <- "group_data"
  
  # Two-sample t-test
  fit <- fmri_ttest(gd, formula = ~ 1 + group, engine = "classic")
  
  expect_s3_class(fit, "fmri_ttest_fit")
  expect_equal(nrow(fit$beta), 2)  # Intercept + group
  expect_true("(Intercept)" %in% rownames(fit$beta))
  expect_true("groupB" %in% rownames(fit$beta))
  
  # Group effect should be detected
  group_p <- fit$p["groupB", ]
  expect_true(mean(group_p < 0.05) > 0.3)  # At least 30% significant
})

test_that("fmri_ttest works with Welch t-test", {
  set.seed(125)
  
  # Simulate data with unequal variances
  S <- 20
  P <- 30
  group <- factor(rep(c("A", "B"), each = 10))
  
  Y <- matrix(0, nrow = S, ncol = P)
  Y[group == "A", ] <- matrix(rnorm(10 * P, mean = 0, sd = 1), nrow = 10)
  Y[group == "B", ] <- matrix(rnorm(10 * P, mean = 0.5, sd = 2), nrow = 10)
  
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = NULL,
        covars = data.frame(subject = 1:S, group = group),
        feature = NULL
      )
    )
  )
  class(gd) <- "group_data"
  
  # Welch t-test
  fit <- fmri_ttest(gd, formula = ~ 1 + group, engine = "welch")
  
  expect_s3_class(fit, "fmri_ttest_fit")
  expect_equal(nrow(fit$beta), 2)
  expect_true(all(is.finite(fit$t["group", ])))
  expect_true(all(fit$df["group", ] < (S - 2)))  # Welch df < pooled df
})

test_that("fmri_ttest works with meta-analysis engine", {
  set.seed(126)
  
  # Simulate data with standard errors
  S <- 15
  P <- 40
  Y <- matrix(rnorm(S * P, mean = 0.3), nrow = S, ncol = P)
  V <- matrix(runif(S * P, 0.1, 0.5), nrow = S, ncol = P)  # Variances
  
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = V,
        covars = data.frame(subject = 1:S),
        feature = NULL
      )
    )
  )
  class(gd) <- "group_data"
  
  # Meta-analysis
  fit <- fmri_ttest(gd, formula = ~ 1, engine = "meta")
  
  expect_s3_class(fit, "fmri_ttest_fit")
  expect_equal(fit$engine, "meta")
  expect_equal(dim(fit$beta), c(1, P))
  expect_equal(dim(fit$se), c(1, P))
  expect_equal(dim(fit$z), c(1, P))
  expect_true(all(fit$df == Inf))  # Meta-analysis uses normal approximation
})

test_that("paired_diff_block works correctly", {
  set.seed(127)
  
  S <- 10
  P <- 20
  
  blkA <- list(
    Y = matrix(rnorm(S * P, mean = 1), nrow = S, ncol = P),
    V = matrix(runif(S * P, 0.1, 0.3), nrow = S, ncol = P),
    meta = list(subjects = paste0("sub", 1:S), contrast = "A"),
    covars = data.frame(subject = 1:S)
  )
  
  blkB <- list(
    Y = matrix(rnorm(S * P, mean = 0.5), nrow = S, ncol = P),
    V = matrix(runif(S * P, 0.1, 0.3), nrow = S, ncol = P),
    meta = list(subjects = paste0("sub", 1:S), contrast = "B"),
    covars = data.frame(subject = 1:S)
  )
  
  # Compute difference
  diff_blk <- paired_diff_block(blkA, blkB, rho = 0)
  
  expect_equal(dim(diff_blk$Y), c(S, P))
  expect_equal(diff_blk$Y, blkA$Y - blkB$Y)
  expect_equal(diff_blk$V, blkA$V + blkB$V)  # rho = 0
  expect_equal(diff_blk$meta$contrast, "A_minus_B")
  
  # With correlation
  diff_blk_corr <- paired_diff_block(blkA, blkB, rho = 0.5)
  expect_true(all(diff_blk_corr$V < diff_blk$V))  # Positive correlation reduces variance
})

test_that("fmri_ttest handles multiple comparisons correction", {
  set.seed(128)
  
  # Simulate data with some true positives
  S <- 20
  P <- 100
  Y <- matrix(rnorm(S * P), nrow = S, ncol = P)
  Y[, 1:10] <- Y[, 1:10] + 2  # Strong signal in first 10 features
  
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = NULL,
        covars = data.frame(subject = 1:S),
        feature = list(group = rep(1:10, each = 10))
      )
    )
  )
  class(gd) <- "group_data"
  
  # Without correction
  fit_nocorr <- fmri_ttest(gd, formula = ~ 1, engine = "classic", mc = NULL)
  
  # With BH correction
  fit_bh <- fmri_ttest(gd, formula = ~ 1, engine = "classic", mc = "bh")
  
  # With spatial FDR
  fit_sfdr <- fmri_ttest(gd, formula = ~ 1, engine = "classic", 
                        mc = "spatial_fdr", alpha = 0.05)
  
  expect_null(fit_nocorr$q)
  expect_equal(dim(fit_bh$q), dim(fit_bh$p))
  expect_equal(dim(fit_sfdr$q), dim(fit_sfdr$p))
  
  # BH q-values should be >= p-values
  expect_true(all(fit_bh$q >= fit_bh$p))
  
  # Should detect true positives
  expect_true(mean(fit_sfdr$q[1, 1:10] < 0.05) > 0.5)
})

test_that("flip_sign works correctly", {
  set.seed(129)
  
  fit <- list(
    beta = matrix(c(1, -2, 3), nrow = 3, 
                  dimnames = list(c("(Intercept)", "groupB", "age"), NULL)),
    t = matrix(c(2, -4, 6), nrow = 3,
              dimnames = list(c("(Intercept)", "groupB", "age"), NULL)),
    z = matrix(c(1.5, -3, 4.5), nrow = 3,
              dimnames = list(c("(Intercept)", "groupB", "age"), NULL))
  )
  
  # Flip all
  fit_flipped <- flip_sign(fit)
  expect_equal(fit_flipped$beta, -fit$beta)
  expect_equal(fit_flipped$t, -fit$t)
  expect_equal(fit_flipped$z, -fit$z)
  
  # Flip specific coefficient
  fit_partial <- flip_sign(fit, coef = "groupB")
  expect_equal(fit_partial$beta["(Intercept)", ], fit$beta["(Intercept)", ])
  expect_equal(fit_partial$beta["groupB", ], -fit$beta["groupB", ])
  expect_equal(fit_partial$beta["age", ], fit$beta["age", ])
})

test_that("print and summary methods work", {
  set.seed(130)
  
  S <- 10
  P <- 20
  Y <- matrix(rnorm(S * P), nrow = S, ncol = P)
  
  gd <- list(
    blocks = list(
      list(
        Y = Y,
        V = NULL,
        covars = data.frame(subject = 1:S),
        feature = NULL
      )
    )
  )
  class(gd) <- "group_data"
  
  fit <- fmri_ttest(gd, formula = ~ 1, engine = "classic")
  
  expect_output(print(fit), "fmri_ttest Results")
  expect_output(print(fit), "Engine: classic")
  expect_output(summary(fit), "Coefficients")
})