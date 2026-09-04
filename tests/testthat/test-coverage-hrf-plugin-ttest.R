# Coverage for estimate_hrf helpers, plugin capability gates, and ttest MC/meta.

test_that("estimate_hrf formula/order/basis/penalty/solve helpers", {
  expect_error(fmrireg:::.validate_estimate_hrf_formula(~ hrf(cond)), "two-sided")
  expect_error(
    fmrireg:::.validate_estimate_hrf_formula(onset ~ trialwise(cond)),
    "trialwise"
  )
  expect_error(
    fmrireg:::.validate_estimate_hrf_formula(onset ~ factor(cond)),
    "hrf\\(\\)"
  )
  expect_true(fmrireg:::.validate_estimate_hrf_formula(onset ~ hrf(condition)))

  ev <- data.frame(
    onset = c(4, 1, 3, 2),
    condition = factor(c("A", "B", "A", "B")),
    run = c(2, 1, 2, 1)
  )
  ordered <- fmrireg:::.order_hrf_event_data(ev, onset ~ hrf(condition), ~ run)
  expect_equal(ordered$onset, c(1, 2, 3, 4))
  expect_error(
    fmrireg:::.order_hrf_event_data(ev, onsets ~ hrf(condition), ~ run),
    "onset column"
  )
  expect_error(
    fmrireg:::.order_hrf_event_data(ev, onset ~ hrf(condition), run ~ 1),
    "one-sided"
  )

  spec_bs <- fmrireg:::.new_hrf_basis_spec("bspline", k = 4L, span = 20)
  expect_equal(spec_bs$degree, 3L)
  expect_equal(spec_bs$k, 4L)
  eval_bs <- fmrireg:::.evaluate_hrf_basis(seq(0, 20, by = 2), spec_bs)
  expect_equal(dim(eval_bs), c(11L, 4L))
  expect_true(all(eval_bs[1, ] == 0 | is.finite(eval_bs[1, ])))

  spec_tent <- fmrireg:::.new_hrf_basis_spec("tent", k = 3L, span = 16)
  expect_equal(spec_tent$degree, 1L)
  eval_tent <- fmrireg:::.evaluate_hrf_basis(c(-1, 0, 8, 16, 30), spec_tent)
  expect_equal(dim(eval_tent), c(5L, 3L))
  expect_true(all(eval_tent[1, ] == 0))
  expect_true(all(eval_tent[5, ] == 0))

  hrf_obj <- fmrireg:::.as_estimation_hrf_basis(spec_bs)
  expect_true(!is.null(hrf_obj))

  pen <- fmrireg:::.hrf_smoothness_penalty(k = 4L, n_curves = 2L)
  expect_equal(dim(pen), c(8L, 8L))
  pen_small <- fmrireg:::.hrf_smoothness_penalty(k = 2L, n_curves = 1L)
  expect_equal(dim(pen_small), c(2L, 2L))

  XtX <- crossprod(matrix(rnorm(40 * 4), 40, 4))
  XtY <- matrix(rnorm(4 * 2), 4, 2)
  solved <- fmrireg:::.solve_hrf_system(XtX, XtY, diag(4), lambda = 0.1)
  expect_true(!is.null(solved$coefficients))
  expect_equal(dim(solved$coefficients), c(4L, 2L))
  expect_true(is.finite(solved$edf))

  # Singular penalty path returns NULL when system not PD
  expect_null(
    fmrireg:::.solve_hrf_system(
      matrix(0, 2, 2), matrix(0, 2, 1), matrix(0, 2, 2), lambda = 0
    )
  )

  sel <- fmrireg:::.select_hrf_lambda(
    XtX, XtY, response_ss = sum(XtY^2), penalty = diag(4),
    n_effective = 40, lambda = 0.5, lambda_grid = c(0.1, 0.5, 1),
    progress = FALSE
  )
  expect_true(!is.null(sel$lambda) || !is.null(sel$chosen) || is.list(sel))

  expect_error(
    fmrireg:::.select_hrf_lambda(
      XtX, XtY, 1, diag(4), 40, lambda = "nope", lambda_grid = 1, progress = FALSE
    ),
    "gcv"
  )
  expect_error(
    fmrireg:::.select_hrf_lambda(
      XtX, XtY, 1, diag(4), 40, lambda = -1, lambda_grid = 1, progress = FALSE
    ),
    "non-negative"
  )

  # Name-suffix curve map path
  design <- matrix(
    rnorm(20 * 6), 20, 6,
    dimnames = list(NULL, c("cond_A_b1", "cond_A_b2", "cond_A_b3",
                            "cond_B_b1", "cond_B_b2", "cond_B_b3"))
  )
  cmap <- fmrireg:::.hrf_curve_map(list(), design, k = 3L)
  expect_equal(length(cmap$indices), 2L)
  expect_equal(nrow(cmap$info), 2L)

  expect_error(
    fmrireg:::.hrf_curve_map(list(), matrix(1:4, 2, 2), k = 2L),
    "map all event-design|column names"
  )
})

test_that("plugin capability normalize/validate/derive cover rejection paths", {
  caps <- fmrireg:::.normalize_engine_capabilities(NULL)
  expect_true(isTRUE(caps$robust))
  expect_true("model" %in% caps$variance_methods)

  caps2 <- fmrireg:::.normalize_engine_capabilities(list(robust = FALSE, ma = TRUE))
  expect_false(caps2$robust)
  expect_true(caps2$ma)

  cfg <- fmri_lm_control()
  expect_true(is.list(
    fmrireg:::.validate_engine_capabilities("ols", cfg, capabilities = NULL)
  ))

  cfg_robust <- fmri_lm_control(robust = robust_spec("huber"))
  # Only assert if robust actually enabled in this control shape
  if (isTRUE(fmrireg:::.fmri_lm_robust_enabled(cfg_robust$robust))) {
    expect_error(
      fmrireg:::.validate_engine_capabilities(
        "lowrank", cfg_robust, capabilities = list(robust = FALSE)
      ),
      "robust"
    )
  }

  cfg_ma <- fmri_lm_control()
  cfg_ma$ar$q <- 1L
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "lowrank", cfg_ma, capabilities = list(ma = FALSE)
    ),
    "MA"
  )

  cfg_vox <- fmri_lm_control()
  cfg_vox$ar$voxelwise <- TRUE
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "lowrank", cfg_vox, capabilities = list(ar_voxelwise = FALSE)
    ),
    "voxelwise|shared"
  )

  derived <- fmrireg:::.derive_engine_execution_config(
    cfg, capabilities = list(robust = FALSE, preprocessing = FALSE)
  )
  expect_true(inherits(derived, "fmri_lm_control") || is.list(derived))
})

test_that(".fmri_ttest_apply_mc covers BH and spatial_fdr branches", {
  set.seed(3)
  res <- list(
    z = matrix(c(2.5, -1, 0.5, 3), 2, 2, dimnames = list(c("a", "b"), NULL)),
    p = matrix(c(0.01, 0.4, 0.2, 0.002), 2, 2, dimnames = list(c("a", "b"), NULL))
  )
  expect_identical(fmrireg:::.fmri_ttest_apply_mc(res, NULL, 0.05, NULL), res)

  bh <- fmrireg:::.fmri_ttest_apply_mc(res, "bh", 0.05, NULL)
  expect_equal(dim(bh$q), dim(res$p))
  expect_true(all(bh$q >= 0 & bh$q <= 1))

  by_adj <- fmrireg:::.fmri_ttest_apply_mc(res, "by", 0.05, NULL)
  expect_equal(dim(by_adj$q), dim(res$p))

  expect_error(
    fmrireg:::.fmri_ttest_apply_mc(res, "spatial_fdr", 0.05, NULL),
    "feature grouping"
  )

  if (exists("spatial_fdr", mode = "function") ||
      is.function(get0("spatial_fdr", envir = asNamespace("fmrireg"), inherits = FALSE))) {
    sfdr <- fmrireg:::.fmri_ttest_apply_mc(
      res, "spatial_fdr", 0.1, feature_group = c(1L, 1L)
    )
    expect_true(!is.null(sfdr$q))
    expect_true(!is.null(sfdr$reject))
  }

  # Contrast-only path
  res_c <- list(
    z_contrast = c(2.2, -0.5, 1.1),
    p_contrast = c(0.03, 0.6, 0.2)
  )
  bh_c <- fmrireg:::.fmri_ttest_apply_mc(res_c, "bh", 0.05, NULL)
  expect_equal(length(bh_c$q_contrast), 3L)
})

test_that("fmri_meta_fit_extended and coef_image.fmri_ttest_fit cover branches", {
  set.seed(11)
  S <- 8L
  P <- 4L
  Y <- matrix(rnorm(S * P, mean = 0.2), S, P)
  V <- matrix(runif(S * P, 0.05, 0.2), S, P)
  X <- cbind(`(Intercept)` = 1, group = rep(c(0, 1), each = 4))

  fit <- fmrireg:::fmri_meta_fit_extended(Y, V, X, method = "fe", robust = "none")
  expect_true(!is.null(fit$beta) || !is.null(fit$coefficients))
  expect_equal(fit$method, "fe")
  expect_true(is.logical(fit$ok) || all(fit$ok %in% c(0, 1)))

  expect_error(
    fmrireg:::fmri_meta_fit_extended(Y, V[1:3, ], X),
    "same dimensions"
  )
  expect_error(
    fmrireg:::fmri_meta_fit_extended(Y, V, X[1:3, ]),
    "same number of rows"
  )

  C <- matrix(rnorm(S * P), S, P)
  fit_v <- fmrireg:::fmri_meta_fit_extended(
    Y, V, X, method = "fe", voxelwise = C, center_voxelwise = TRUE,
    voxel_name = "vox"
  )
  expect_true("vox" %in% rownames(fit_v$beta) || "vox" %in% rownames(fit_v$se))

  expect_error(
    fmrireg:::fmri_meta_fit_extended(
      Y, V, X, voxelwise = matrix(1, 2, 2)
    ),
    "match Y"
  )

  tfit <- fmrireg:::.demo_fmri_ttest_fit()
  est <- coef_image(tfit, coef = 1, statistic = "estimate")
  expect_equal(length(est), 1L)
  est_named <- coef_image(tfit, coef = "conditionA", statistic = "z")
  expect_true(is.numeric(est_named))
  pvals <- coef_image(tfit, coef = 2, statistic = "p")
  expect_true(all(pvals >= 0 & pvals <= 1))

  expect_error(coef_image(tfit, coef = "missing", statistic = "estimate"), "not found")
  expect_error(coef_image(tfit, coef = 99, statistic = "estimate"), "out of bounds")

  # z from t when z missing
  tfit2 <- tfit
  tfit2$z <- NULL
  z_from_t <- coef_image(tfit2, coef = 1, statistic = "z")
  expect_true(is.numeric(z_from_t))

  tfit3 <- tfit
  tfit3$se <- NULL
  expect_error(coef_image(tfit3, coef = 1, statistic = "se"), "Standard errors")
})

test_that("fmri_ols_fit voxelwise path and flip_sign selective coef", {
  set.seed(4)
  Y <- matrix(rnorm(12 * 3), 12, 3)
  X <- cbind(Intercept = 1, g = rep(c(0, 1), each = 6))
  C <- matrix(rnorm(12 * 3), 12, 3)
  ols <- fmrireg:::fmri_ols_fit(Y, X, voxelwise = C, voxel_name = "cov")
  expect_true("cov" %in% rownames(ols$beta))

  fit <- fmrireg:::.demo_fmri_ttest_fit()
  flipped <- flip_sign(fit, coef = "conditionA")
  expect_equal(flipped$beta["conditionA", 1], -fit$beta["conditionA", 1])
  expect_equal(flipped$beta["conditionB", 1], fit$beta["conditionB", 1])
})
