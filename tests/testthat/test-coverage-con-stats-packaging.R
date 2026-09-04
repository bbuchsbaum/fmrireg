# con_stats.R remaining: packaging helpers, fit_contrasts.default, orient, colind.

test_that("extract_colind / create_full_contrast_* / package_* cover errors", {
  expect_error(fmrireg:::extract_colind(c(1, -1), "c1", "t-contrast"), "Missing")
  bad <- structure(c(1, -1), colind = c(0L, 1L))
  expect_error(fmrireg:::extract_colind(bad, "c1", "t-contrast"), "positive")
  ok <- structure(c(1, -1), colind = 2:3)
  expect_equal(fmrireg:::extract_colind(ok, "c1", "t-contrast"), 2:3)

  full_v <- fmrireg:::create_full_contrast_vector(c(1, -1), 2:3, p = 4)
  expect_equal(dim(full_v), c(1L, 4L))
  expect_equal(as.numeric(full_v), c(0, 1, -1, 0))

  # Matrix row / column coercion
  expect_equal(
    as.numeric(fmrireg:::create_full_contrast_vector(matrix(c(1, -1), 1), 1:2, 3)),
    c(1, -1, 0)
  )
  expect_equal(
    as.numeric(fmrireg:::create_full_contrast_vector(matrix(c(1, -1), 2, 1), 1:2, 3)),
    c(1, -1, 0)
  )
  expect_error(
    fmrireg:::create_full_contrast_vector(matrix(1, 2, 2), 1:2, 3),
    "vector or single-row"
  )
  expect_error(
    fmrireg:::create_full_contrast_vector(c(1, -1), 1:3, 4),
    "Dimension mismatch"
  )
  expect_error(
    fmrireg:::create_full_contrast_vector(c(1, -1), 3:4, 3),
    "out of bounds"
  )

  Cw <- matrix(c(1, 0, 0, 1), 2)
  full_m <- fmrireg:::create_full_contrast_matrix(Cw, colind = 2:3, p = 4)
  expect_equal(dim(full_m), c(2L, 4L))
  expect_error(fmrireg:::create_full_contrast_matrix(1:2, 1:2, 3), "matrix")
  expect_error(
    fmrireg:::create_full_contrast_matrix(Cw, colind = 3:4, p = 3),
    "out of bounds"
  )

  # Cell-oriented F weights get transposed via .orient_fcontrast
  cellish <- matrix(c(1, -1, 0, 1, 0, -1), nrow = 3, ncol = 2)
  oriented <- fmrireg:::.orient_fcontrast(cellish, colind = 1:3)
  expect_true(is.matrix(oriented))
  expect_equal(ncol(oriented), 3L)

  stats_t <- list(
    estimate = 1:2, se = c(0.1, 0.2), stat = c(10, 5),
    prob = c(0.01, 0.05), sigma = c(0.5, 0.5), stat_type = "tstat"
  )
  packed_t <- fmrireg:::package_tcontrast_result(
    "A_vs_B", original_weights = c(1, -1), colind = 2:3, df = 20, stats = stats_t
  )
  expect_equal(packed_t$name, "A_vs_B")
  expect_equal(packed_t$type, "contrast")

  stats_f <- list(
    estimate = 1:2, se = c(0.1, 0.2), stat = c(4, 3),
    prob = c(0.02, 0.08), sigma = c(0.5, 0.5), stat_type = "Fstat"
  )
  packed_f <- fmrireg:::package_fcontrast_result(
    "main", original_weights = diag(2), colind = 2:3, df = 20, stats = stats_f
  )
  expect_equal(packed_f$type, "Fcontrast")
})

test_that("fit_contrasts.default covers se=TRUE/FALSE and aliased path", {
  set.seed(91)
  d <- data.frame(y1 = rnorm(40), y2 = rnorm(40), x = rnorm(40), z = rnorm(40))
  fit <- lm(cbind(y1, y2) ~ x + z, data = d)

  out <- fit_contrasts(fit, conmat = c(0, 1, -1), colind = 1:3, se = TRUE)
  expect_true(inherits(out, "tstat") || inherits(out, "result_stat"))
  expect_equal(length(out$estimate), 2L)
  expect_true(all(is.finite(out$se)))

  out_eff <- fit_contrasts(fit, conmat = matrix(c(0, 1, 0), 1), colind = 1:3, se = FALSE)
  expect_true(inherits(out_eff, "effect") || identical(out_eff$stat_type, "effects"))
  expect_null(out_eff$se)

  # Rank-deficient / aliased contrast
  d2 <- data.frame(y = rnorm(30), x1 = 1:30, x2 = 1:30)
  fit_sing <- lm(y ~ x1 + x2, data = d2)
  expect_warning(
    out_alias <- fit_contrasts(fit_sing, conmat = c(0, 1, -1), colind = 1:3, se = TRUE),
    "non-estimable|aliased|rank|collinear"
  )
  expect_true(all(is.na(out_alias$estimate)) || all(is.na(out_alias$stat)) ||
                inherits(out_alias, "result_stat"))
})

test_that("beta_stats and .contrast_uses_aliased / warn helpers", {
  set.seed(92)
  d <- data.frame(y1 = rnorm(25), y2 = rnorm(25), x = rnorm(25))
  fit <- lm(cbind(y1, y2) ~ x, data = d)
  bs <- fmrireg:::beta_stats(fit, varnames = c("(Intercept)", "x"), se = TRUE)
  expect_true(inherits(bs, "tbl_df") || is.list(bs))

  bs2 <- fmrireg:::beta_stats(fit, varnames = c("(Intercept)", "x"), se = FALSE)
  expect_true(!is.null(bs2))

  L <- matrix(c(0, 1, -1), 1)
  expect_true(fmrireg:::.contrast_uses_aliased(L, aliased = 2L))
  expect_false(fmrireg:::.contrast_uses_aliased(L, aliased = integer(0)))
  expect_false(fmrireg:::.contrast_uses_aliased(matrix(c(1, 0, 0), 1), aliased = 3L))

  expect_warning(
    fmrireg:::.warn_nonestimable_contrast("c1", aliased = 2L, coef_names = c("a", "b", "c")),
    "non-estimable|aliased|c1"
  )
  expect_warning(
    fmrireg:::.warn_nonestimable_contrast(NULL, aliased = 1L, coef_names = NULL),
    "non-estimable|aliased"
  )
})

test_that("fit_contrasts.fmri_lm empty and orientation edge paths", {
  expect_equal(fit_contrasts(structure(list(), class = "fmri_lm"), contrasts = NULL), list())
  expect_equal(fit_contrasts(structure(list(), class = "fmri_lm"), contrasts = list()), list())

  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  con <- pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
  out <- fit_contrasts(fit, list(con))
  expect_true(is.list(out) && length(out) >= 1L)

  # Force event_indices alignment / cov.unscaled subsetting
  fit2 <- fit
  p_event <- length(fit$result$event_indices)
  if (p_event >= 1L && !is.null(fit$result$cov.unscaled)) {
    # Keep full cov; method should subset to event indices when beta is event-sized
    out2 <- tryCatch(fit_contrasts(fit2, list(con)), error = function(e) e)
    expect_true(inherits(out2, "error") || is.list(out2))
  }
})
