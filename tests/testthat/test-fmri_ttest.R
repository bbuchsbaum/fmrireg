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
  expect_true("group" %in% rownames(fit$beta))
  
  # Group effect should be detected
  group_p <- fit$p["group", ]
  expect_true(mean(group_p < 0.05) > 0.3)  # At least 30% significant
  expect_true(mean(fit$beta["group", ]) < 0) # default sign is A - B
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

  fit$p_contrast <- c(0.1, 0.2)
  fit$t_contrast <- c(-2, 3)
  fit_contrast <- flip_sign(fit)
  expect_equal(fit_contrast$p_contrast, fit$p_contrast)
  expect_equal(fit_contrast$t_contrast, -fit$t_contrast)
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

test_that("group coefficient sign convention is stable across engines", {
  set.seed(131)

  S <- 24
  P <- 10
  group <- factor(rep(c("A", "B"), each = S / 2))
  Y <- matrix(rnorm(S * P, sd = 0.1), nrow = S, ncol = P)
  Y[group == "B", ] <- Y[group == "B", ] + 2
  V <- matrix(0.2, nrow = S, ncol = P)

  gd <- structure(
    list(
      blocks = list(
        list(
          Y = Y,
          V = V,
          covars = data.frame(group = group),
          feature = NULL
        )
      )
    ),
    class = "group_data"
  )

  fit_classic <- fmri_ttest(gd, formula = ~ 1 + group, engine = "classic")
  fit_classic_flip <- fmri_ttest(gd, formula = ~ 1 + group, engine = "classic", sign = "BminusA")
  fit_welch <- fmri_ttest(gd, formula = ~ 1 + group, engine = "welch")
  fit_welch_flip <- fmri_ttest(gd, formula = ~ 1 + group, engine = "welch", sign = "BminusA")
  fit_meta <- fmri_ttest(gd, formula = ~ 1 + group, engine = "meta")
  fit_meta_flip <- fmri_ttest(gd, formula = ~ 1 + group, engine = "meta", sign = "BminusA")

  expect_equal(rownames(fit_classic$beta), c("(Intercept)", "group"))
  expect_equal(rownames(fit_welch$beta), c("(Intercept)", "group"))
  expect_equal(rownames(fit_meta$beta), c("(Intercept)", "group"))
  expect_true(mean(fit_classic$beta["group", ]) < 0)
  expect_true(mean(fit_welch$beta["group", ]) < 0)
  expect_true(mean(fit_meta$beta["group", ]) < 0)
  expect_equal(fit_classic_flip$beta["group", ], -fit_classic$beta["group", ])
  expect_equal(fit_welch_flip$beta["group", ], -fit_welch$beta["group", ])
  expect_equal(fit_meta_flip$beta["group", ], -fit_meta$beta["group", ])
})

test_that("fmri_ttest computes exact contrasts where supported", {
  set.seed(132)

  S <- 18
  P <- 8
  group <- factor(rep(c("A", "B"), each = S / 2))
  age <- seq(20, 54, length.out = S)
  Y <- matrix(rnorm(S * P, sd = 0.15), nrow = S, ncol = P)
  Y[group == "B", ] <- Y[group == "B", ] + 1
  Y <- Y + outer(as.numeric(scale(age)), rep(0.1, P))
  V <- matrix(0.3, nrow = S, ncol = P)

  gd <- structure(
    list(
      blocks = list(
        list(
          Y = Y,
          V = V,
          covars = data.frame(group = group, age = age),
          feature = NULL
        )
      )
    ),
    class = "group_data"
  )

  contrast <- c(group = 1, age = 0.5)
  fit_classic <- fmri_ttest(gd, formula = ~ 1 + group + age, engine = "classic", contrast = contrast)
  fit_meta <- fmri_ttest(gd, formula = ~ 1 + group + age, engine = "meta", contrast = contrast)

  expect_true(all(is.finite(fit_classic$z_contrast)))
  expect_true(all(is.finite(fit_meta$z_contrast)))
  expect_equal(length(fit_classic$z_contrast), P)
  expect_equal(length(fit_meta$z_contrast), P)
  expect_true(all(fit_classic$p_contrast >= 0 & fit_classic$p_contrast <= 1))
  expect_true(all(fit_meta$p_contrast >= 0 & fit_meta$p_contrast <= 1))
})

test_that("fmri_ttest rejects unsupported contrast combinations explicitly", {
  set.seed(133)

  S <- 12
  P <- 6
  group <- factor(rep(c("A", "B"), each = S / 2))
  Y <- matrix(rnorm(S * P), nrow = S, ncol = P)
  V <- matrix(0.2, nrow = S, ncol = P)
  voxel_cov <- matrix(rnorm(S * P), nrow = S, ncol = P)

  gd <- structure(
    list(
      blocks = list(
        list(
          Y = Y,
          V = V,
          covars = data.frame(group = group),
          feature = NULL
        )
      )
    ),
    class = "group_data"
  )

  expect_error(
    fmri_ttest(gd, formula = ~ 1 + group, engine = "meta", contrast = c(group = 1), voxelwise_cov = voxel_cov),
    "contrast is not supported with voxelwise_cov in meta analyses"
  )
  expect_error(
    fmri_ttest(gd, formula = ~ 1 + group, engine = "welch", contrast = c(`(Intercept)` = 1, group = 1)),
    "Welch contrasts currently support only the group coefficient"
  )
})

test_that("fmri_ttest errors when spatial_fdr metadata is unavailable", {
  set.seed(134)

  S <- 15
  P <- 20
  gd <- structure(
    list(
      blocks = list(
        list(
          Y = matrix(rnorm(S * P), nrow = S, ncol = P),
          V = NULL,
          covars = data.frame(subject = seq_len(S)),
          feature = NULL
        )
      )
    ),
    class = "group_data"
  )

  expect_error(
    fmri_ttest(gd, formula = ~ 1, engine = "classic", mc = "spatial_fdr"),
    "requires feature grouping metadata"
  )
})

test_that("fmri_ttest validates custom weights", {
  set.seed(135)

  S <- 10
  P <- 5
  gd <- structure(
    list(
      blocks = list(
        list(
          Y = matrix(rnorm(S * P), nrow = S, ncol = P),
          V = matrix(0.2, nrow = S, ncol = P),
          covars = data.frame(subject = seq_len(S)),
          feature = NULL
        )
      )
    ),
    class = "group_data"
  )

  expect_error(
    fmri_ttest(gd, formula = ~ 1, engine = "meta", weights = "custom", weights_custom = rep(1, S - 1)),
    "weights_custom must be length S or SxP"
  )
  expect_error(
    fmri_ttest(gd, formula = ~ 1, engine = "meta", weights = "custom", weights_custom = c(rep(1, S - 1), 0)),
    "finite positive"
  )
})

test_that("fmri_ttest covers paired, auto, mu0, voxelwise_cov, and coef_image", {
  set.seed(136)

  S <- 14
  P <- 12
  group <- factor(rep(c("A", "B"), each = S / 2))
  Y <- matrix(rnorm(S * P, mean = 2), nrow = S, ncol = P)
  V <- matrix(0.25, nrow = S, ncol = P)
  voxel_cov <- matrix(rnorm(S * P), nrow = S, ncol = P)

  gd <- structure(
    list(
      blocks = list(
        list(
          Y = Y,
          V = V,
          covars = data.frame(group = group),
          feature = list(group = rep(1:3, each = 4))
        )
      )
    ),
    class = "group_data"
  )

  fit_auto <- fmri_ttest(gd, formula = ~ 1, engine = "auto")
  fit_mu0 <- fmri_ttest(gd, formula = ~ 1, engine = "classic", mu0 = 2)
  fit_paired <- fmri_ttest(gd, formula = ~ 1, engine = "classic", paired = TRUE)
  fit_vcov <- fmri_ttest(gd, formula = ~ 1, engine = "classic", voxelwise_cov = voxel_cov)

  expect_equal(fit_auto$engine, "meta")
  expect_true(abs(mean(fit_mu0$beta["(Intercept)", ])) < 0.5)
  expect_equal(rownames(fit_paired$beta), "(Intercept)")
  expect_true("voxel_cov" %in% rownames(fit_vcov$beta))
  expect_equal(coef_image(fit_vcov, coef = "voxel_cov", statistic = "estimate"), fit_vcov$beta["voxel_cov", ])
  expect_true(all(is.finite(coef_image(fit_vcov, coef = "voxel_cov", statistic = "z"))))
})
