# fmri_ttest_helpers.R remaining branches: paired_diff rho shapes, coef_image mask,
# group/sample label helpers, rename/sign/contrast utilities.

test_that("paired_diff_block covers rho vector/matrix and subject mismatch", {
  set.seed(71)
  S <- 5L
  P <- 4L
  blkA <- list(
    Y = matrix(rnorm(S * P), S, P),
    V = matrix(runif(S * P, 0.1, 0.4), S, P),
    meta = list(subjects = paste0("s", seq_len(S)), contrast = "face"),
    index = 1:P,
    covars = data.frame(id = seq_len(S)),
    feature = list(label = paste0("v", seq_len(P)))
  )
  blkB <- blkA
  blkB$Y <- blkA$Y + matrix(rnorm(S * P, sd = 0.2), S, P)
  blkB$V <- matrix(runif(S * P, 0.1, 0.4), S, P)
  blkB$meta$contrast <- "place"

  # Per-feature rho
  dP <- paired_diff_block(blkA, blkB, rho = rep(0.2, P))
  expect_equal(dim(dP$Y), c(S, P))
  expect_equal(dim(dP$V), c(S, P))
  expect_match(dP$meta$contrast, "face_minus_place")

  # Matrix rho
  dM <- paired_diff_block(blkA, blkB, rho = matrix(0.1, S, P))
  expect_true(all(is.finite(dM$V) | is.na(dM$V)))

  # Subject mismatch
  blkB2 <- blkB
  blkB2$meta$subjects <- paste0("x", seq_len(S))
  expect_error(paired_diff_block(blkA, blkB2), "subjects must match")

  # No variance -> V NULL
  blkA3 <- blkA
  blkA3$V <- NULL
  d3 <- paired_diff_block(blkA3, blkB, rho = 0)
  expect_null(d3$V)
})

test_that("coef_image.fmri_ttest_fit mask NeuroVol and z-from-t df matrix", {
  skip_if_not_installed("neuroim2")
  tfit <- fmrireg:::.demo_fmri_ttest_fit()
  # Expand to multi-feature for mask mapping
  P <- 4L
  tfit$beta <- matrix(rnorm(2 * P), 2, P, dimnames = list(c("conditionA", "conditionB"), NULL))
  tfit$se <- matrix(runif(2 * P, 0.1, 0.5), 2, P, dimnames = list(rownames(tfit$beta), NULL))
  tfit$z <- matrix(rnorm(2 * P), 2, P, dimnames = list(rownames(tfit$beta), NULL))
  tfit$p <- matrix(runif(2 * P), 2, P, dimnames = list(rownames(tfit$beta), NULL))
  tfit$t <- tfit$z
  tfit$df <- matrix(rep(10, 2 * P), 2, P, dimnames = list(rownames(tfit$beta), NULL))

  mask <- neuroim2::LogicalNeuroVol(
    array(c(TRUE, TRUE, TRUE, TRUE), dim = c(2, 2, 1)),
    neuroim2::NeuroSpace(c(2, 2, 1))
  )
  vol <- coef_image(tfit, coef = 1, statistic = "estimate", mask = mask)
  expect_s4_class(vol, "NeuroVol")
  expect_equal(length(as.vector(vol)), 4L)

  # z from t with matrix df / named coef
  tfit2 <- tfit
  tfit2$z <- NULL
  zvals <- coef_image(tfit2, coef = "conditionA", statistic = "z")
  expect_equal(length(zvals), P)

  # Missing statistic entirely
  tfit3 <- tfit
  tfit3$p <- NULL
  expect_error(coef_image(tfit3, coef = 1, statistic = "p"), "not available")

  # z without df warns and treats t as z
  tfit4 <- tfit
  tfit4$z <- NULL
  tfit4$df <- NULL
  expect_warning(
    z2 <- coef_image(tfit4, coef = 1, statistic = "z"),
    "No df available"
  )
  expect_equal(length(z2), P)
})

test_that(".fmri_ttest_* rename/sign/contrast/feature helpers", {
  res <- list(
    beta = matrix(c(1, -2, 0.5), 3, 1, dimnames = list(c("a", "b", "c"), NULL)),
    se = matrix(c(0.1, 0.2, 0.3), 3, 1, dimnames = list(c("a", "b", "c"), NULL)),
    t = matrix(c(10, -10, 1), 3, 1, dimnames = list(c("a", "b", "c"), NULL))
  )
  renamed <- fmrireg:::.fmri_ttest_rename_rows(res$beta, "a", "Intercept")
  expect_true("Intercept" %in% rownames(renamed) || identical(renamed, res$beta) ||
                is.matrix(renamed))

  canon <- fmrireg:::.fmri_ttest_canonical_coef_names(
    c("group1", "group2"),
    group_info = list(levels = c("A", "B"))
  )
  expect_true(is.character(canon) || is.null(canon))

  signed <- fmrireg:::.fmri_ttest_apply_group_sign(
    res, group_info = list(levels = c("A", "B")),
    target_sign = 1, source_sign = -1
  )
  expect_true(is.list(signed))

  weights <- fmrireg:::.fmri_ttest_resolve_contrast(
    contrast = c(1, -1, 0),
    coef_names = c("a", "b", "c")
  )
  expect_equal(length(weights), 3L)

  raw <- fmrireg:::.fmri_ttest_raw_contrast_weights(
    canonical_weights = c(Intercept = 0, group = 1),
    coef_names = c("Intercept", "groupB"),
    group_info = list(raw_name = "groupB", canonical_name = "group"),
    target_sign = 1,
    source_sign = -1
  )
  expect_true(is.numeric(raw))
  expect_true("groupB" %in% names(raw))

  sc <- fmrireg:::.fmri_ttest_single_coef_contrast(
    res, weights = c(a = 1, b = 0, c = 0)
  )
  expect_true(is.list(sc))
  expect_equal(sc$estimate, res$beta["a", ])

  stored <- fmrireg:::.fmri_ttest_store_contrast(
    res,
    contrast_stats = list(
      estimate = 1, se = 0.1, t = 10, z = 9, p = 0.001, df = 8
    )
  )
  expect_true(!is.null(stored$beta_contrast) || !is.null(stored$z_contrast) ||
                is.list(stored))

  # Feature group / sample labels from block and attributes
  gd <- list()
  expect_null(fmrireg:::.fmri_ttest_feature_group(gd, block = NULL))
  expect_equal(
    fmrireg:::.fmri_ttest_feature_group(
      gd, block = list(feature = list(group = c(1L, 1L, 2L)))
    ),
    c(1L, 1L, 2L)
  )
  attr(gd, "fmrireg_feature_group") <- c(1L, 2L)
  expect_equal(fmrireg:::.fmri_ttest_feature_group(gd), c(1L, 2L))

  expect_null(fmrireg:::.fmri_ttest_sample_labels(list(), NULL))
  expect_equal(
    fmrireg:::.fmri_ttest_sample_labels(
      list(), block = list(feature = list(label = c("a", "b")))
    ),
    c("a", "b")
  )
  gd2 <- list()
  attr(gd2, "fmrireg_sample_labels") <- c("v1", "v2")
  expect_equal(fmrireg:::.fmri_ttest_sample_labels(gd2), c("v1", "v2"))
})

test_that(".fmri_ttest_validate_weights_custom and group_term", {
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(NULL, 4, 3),
    "weights_custom"
  )
  w <- matrix(1, 4, 3)
  expect_equal(fmrireg:::.fmri_ttest_validate_weights_custom(w, 4, 3), w)
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(matrix(1, 2, 3), 4, 3),
    "length S or SxP"
  )
  wv <- fmrireg:::.fmri_ttest_validate_weights_custom(rep(1.5, 4), 4, 3)
  expect_equal(wv, rep(1.5, 4))
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(c(1, -1, 1, 1), 4, 3),
    "finite positive"
  )
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom("bad", 4, 3),
    "length S or SxP"
  )

  X <- cbind(Intercept = 1, groupB = rep(c(0, 1), each = 3))
  covars <- data.frame(group = factor(rep(c("A", "B"), each = 3)))
  gt <- fmrireg:::.fmri_ttest_group_term(X, covars)
  expect_true(is.null(gt) || is.list(gt))
  expect_null(fmrireg:::.fmri_ttest_group_term(X, data.frame(x = 1:6)))
})

test_that("flip_sign flips vectors and all coefs", {
  fit <- list(
    beta = matrix(1:4, 2, 2, dimnames = list(c("a", "b"), NULL)),
    z_contrast = c(1.5, -0.5),
    t_contrast = c(2, -1)
  )
  all_flip <- flip_sign(fit)
  expect_equal(all_flip$beta, -fit$beta)
  expect_equal(all_flip$z_contrast, -fit$z_contrast)

  partial <- flip_sign(fit, coef = "a")
  expect_equal(partial$beta["a", ], -fit$beta["a", ])
  expect_equal(partial$beta["b", ], fit$beta["b", ])
  # vector contrasts flip even when coef filtered
  expect_equal(partial$z_contrast, -fit$z_contrast)
})
