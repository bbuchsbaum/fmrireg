# Coverage for con_stats.R and fmri_meta / fmri_meta_methods helpers.

tiny_mlm <- function(n = 50L, seed = 5) {
  set.seed(seed)
  data.frame(
    y1 = rnorm(n),
    y2 = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  ) |>
    (\(d) lm(cbind(y1, y2) ~ x1 + x2, data = d))()
}

toy_meta <- function(n_feat = 8L, with_het = TRUE, roi = FALSE) {
  set.seed(8)
  coefs <- cbind(
    `(Intercept)` = rnorm(n_feat, 0.2, 0.1),
    groupB = rnorm(n_feat, 0.5, 0.15)
  )
  se <- cbind(
    `(Intercept)` = runif(n_feat, 0.05, 0.12),
    groupB = runif(n_feat, 0.06, 0.15)
  )
  obj <- list(
    coefficients = coefs,
    se = se,
    method = "DL",
    robust = "none",
    formula = ~ group,
    n_subjects = 12L,
    tau2 = if (with_het) runif(n_feat, 0, 0.05) else NULL,
    I2 = if (with_het) runif(n_feat, 0, 40) else NULL,
    Q = if (with_het) runif(n_feat, 5, 20) else NULL,
    Q_df = if (with_het) rep(10, n_feat) else NULL
  )
  if (roi) {
    obj$n_rois <- n_feat
    obj$roi_names <- paste0("roi", seq_len(n_feat))
    class(obj) <- c("fmri_meta_roi", "fmri_meta")
  } else {
    obj$n_voxels <- n_feat
    class(obj) <- "fmri_meta"
  }
  obj
}

test_that("fit_contrasts.default se=FALSE and estimate_contrast helpers", {
  fit <- tiny_mlm()

  with_se <- fit_contrasts(fit, c(1, -1), colind = 2:3, se = TRUE)
  expect_s3_class(with_se, "tstat")
  expect_equal(length(with_se$estimate), 2L)
  expect_equal(length(with_se$se), 2L)
  expect_true(all(is.finite(with_se$stat)))

  no_se <- fit_contrasts(fit, c(1, -1), colind = 2:3, se = FALSE)
  expect_s3_class(no_se, "effect")
  expect_equal(no_se$stat_type, "effects")
  expect_equal(length(no_se$estimate), 2L)
  expect_null(no_se$se)
  expect_equal(no_se$estimate, with_se$estimate, tolerance = 1e-12)

  con <- structure(list(name = "x1_vs_x2", weights = matrix(c(1, -1), 1)), class = "contrast")
  est <- fmrireg:::estimate_contrast.contrast(con, fit, colind = 2:3)
  expect_s3_class(est, "tbl_df")
  expect_equal(est$type, "contrast")
  expect_equal(est$name, "x1_vs_x2")
  expect_equal(length(est$data[[1]]$estimate), 2L)

  fcon <- structure(list(name = "both", weights = diag(2)), class = "Fcontrast")
  fest <- fmrireg:::estimate_contrast.Fcontrast(fcon, fit, colind = 2:3)
  expect_equal(fest$type, "Fcontrast")
  expect_equal(fest$stat_type, "Fstat")
  expect_true(all(fest$data[[1]]$stat > 0))
})

test_that("contrast scope enforcement and packaging helpers", {
  expect_true(
    fmrireg:::.enforce_contrast_scope(1:2, list(allowed_colind = 1:4), "ok")
  )
  expect_true(fmrireg:::.enforce_contrast_scope(1:2, NULL, "ok"))
  expect_error(
    fmrireg:::.enforce_contrast_scope(
      c(1L, 5L),
      list(allowed_colind = 1:3, reason = "task-only contrasts"),
      "bad_con"
    ),
    "unsupported coefficient indices"
  )
  expect_error(
    fmrireg:::.enforce_contrast_scope(
      9L,
      list(allowed_colind = 1:2),
      "bare"
    ),
    "bare"
  )

  full_t <- fmrireg:::create_full_contrast_vector(c(1, -1), colind = c(2L, 3L), p = 4L)
  expect_equal(dim(full_t), c(1L, 4L))
  expect_equal(as.numeric(full_t), c(0, 1, -1, 0))
  expect_error(
    fmrireg:::create_full_contrast_vector(c(1, -1), colind = 1:3, p = 4L),
    "Dimension mismatch"
  )
  expect_error(
    fmrireg:::create_full_contrast_vector(1, colind = 5L, p = 3L),
    "out of bounds"
  )

  full_f <- fmrireg:::create_full_contrast_matrix(diag(2), colind = 2:3, p = 4L)
  expect_equal(ncol(full_f), 4L)
  expect_error(
    fmrireg:::create_full_contrast_matrix(1:2, colind = 1:2, p = 3L),
    "must be a matrix"
  )

  expect_equal(
    fmrireg:::extract_colind(structure(1, colind = 2:3), "c1", "t-contrast"),
    2:3
  )
  expect_error(
    fmrireg:::extract_colind(1, "c1", "t-contrast"),
    "Missing 'colind'"
  )
  expect_error(
    fmrireg:::extract_colind(structure(1, colind = c(0, 1)), "c1", "t-contrast"),
    "positive integers"
  )

  packed_t <- fmrireg:::package_tcontrast_result(
    "c1", c(1, -1), 1:2, df = 30,
    stats = list(estimate = 1:3, se = rep(0.2, 3), stat = 2:4, prob = rep(0.05, 3))
  )
  expect_equal(packed_t$name, "c1")
  expect_equal(packed_t$type, "contrast")

  packed_f <- fmrireg:::package_fcontrast_result(
    "f1", diag(2), 1:2, df = 30,
    stats = list(estimate = 1:3, se = rep(1, 3), stat = 2:4, prob = rep(0.1, 3))
  )
  expect_equal(packed_f$type, "Fcontrast")
})

test_that("aliased-contrast helpers and beta_stats se=FALSE", {
  expect_false(fmrireg:::.contrast_uses_aliased(matrix(c(1, -1, 0), 1), integer(0)))
  expect_true(fmrireg:::.contrast_uses_aliased(matrix(c(1, 0, 1), 1), aliased = 3L))
  expect_false(fmrireg:::.contrast_uses_aliased(matrix(c(1, -1, 0), 1), aliased = 3L))

  expect_warning(
    fmrireg:::.warn_nonestimable_contrast("c1", aliased = 3L, coef_names = c("a", "b", "c")),
    "non-estimable|aliased|c1"
  )

  fit <- tiny_mlm()
  bs <- fmrireg:::beta_stats(fit, varnames = c("(Intercept)", "x1", "x2"), se = FALSE)
  expect_equal(bs$stat_type, "effect")
  expect_true(is.null(bs$data[[1]]$se[[1]]) || identical(bs$data[[1]]$se[[1]], NULL))
})

test_that("fmri_meta methods cover formula contrasts, pvalues, print/summary, coef_image", {
  meta <- toy_meta()

  expect_equal(dim(coef(meta)), c(8L, 2L))
  expect_equal(dim(se(meta)), c(8L, 2L))
  expect_equal(dim(zscores(meta)), c(8L, 2L))

  p_two <- pvalues(meta, two_tailed = TRUE)
  p_one <- pvalues(meta, two_tailed = FALSE)
  expect_equal(dim(p_one), dim(p_two))
  # For positive z, one-tailed p is half of two-tailed
  z <- zscores(meta)
  pos <- which(z > 0)
  expect_equal(p_one[pos], p_two[pos] / 2, tolerance = 1e-12)

  num_con <- contrast(meta, c(0, 1))
  expect_s3_class(num_con, "fmri_meta_contrast")
  expect_equal(length(num_con$estimate), 8L)

  # parse_contrast_formula handles simple "~ a - b" using plain coef names
  meta_grp <- meta
  colnames(meta_grp$coefficients) <- c("groupold", "groupyoung")
  colnames(meta_grp$se) <- c("groupold", "groupyoung")
  form_con <- contrast(meta_grp, ~ groupold - groupyoung)
  expect_equal(unname(form_con$weights), c(1, -1))

  named_con <- contrast(meta, c(groupB = 1, `(Intercept)` = -1))
  expect_equal(unname(named_con$weights), c(-1, 1))

  expect_error(contrast(meta, c(1, 2, 3)), "length must match")
  expect_error(contrast(meta, "bad"), "Invalid contrast")
  expect_error(
    contrast(meta, c(missing_coef = 1)),
    "Coefficient\\(s\\) not found"
  )

  # Packed-tri exact SE path
  K <- 2L
  tsize <- K * (K + 1) / 2
  tri <- matrix(c(0.01, 0.002, 0.02), nrow = tsize, ncol = 8)
  meta_cov <- meta
  meta_cov$cov <- list(type = "tri", tri = tri)
  con_exact <- contrast(meta_cov, c(1, -1))
  expect_true(all(is.finite(con_exact$se)))

  out <- capture.output(print(meta))
  expect_true(any(grepl("fMRI Meta-Analysis Results", out)))
  expect_true(any(grepl("Heterogeneity", out)))

  meta_vox <- toy_meta(with_het = FALSE)
  out2 <- capture.output(print(meta_vox))
  expect_true(any(grepl("Voxels analyzed", out2)))

  sout <- capture.output(summary(meta, threshold = 0.1))
  expect_true(any(grepl("Coefficients", sout)))
  expect_true(any(grepl("Significant", sout)))

  img_est <- coef_image(meta, coef = "groupB", statistic = "estimate")
  expect_equal(length(img_est), 8L)
  img_p <- coef_image(meta, coef = 1, statistic = "p")
  expect_true(all(img_p >= 0 & img_p <= 1))
  expect_error(coef_image(meta, coef = "nope"), "not found")

  # reconstruct_image for synthetic nifti-backed meta
  meta_nii <- meta
  meta_nii$data <- structure(
    list(dim = c(2L, 2L, 2L), voxel_size = c(2, 2, 2)),
    class = c("group_data_nifti", "group_data")
  )
  meta_nii$voxel_indices <- seq_len(8)
  vol <- fmrireg:::reconstruct_image(as.numeric(meta$coefficients[, 1]), meta_nii)
  expect_true(inherits(vol, "NeuroVol") || is.array(as.array(vol)))

  meta_h5 <- meta
  meta_h5$data <- structure(
    list(
      dim = c(2L, 2L, 2L, 1L),
      mask = array(TRUE, dim = c(2, 2, 2)),
      space = neuroim2::NeuroSpace(c(2, 2, 2))
    ),
    class = c("group_data_h5", "group_data")
  )
  vol_h5 <- fmrireg:::reconstruct_image(as.numeric(meta$coefficients[, 2]), meta_h5)
  expect_true(inherits(vol_h5, "NeuroVol"))
})

test_that("fmri_meta validation / selector / tidy helpers", {
  meta <- toy_meta(n_feat = 5L, roi = TRUE)

  expect_true(fmrireg:::.validate_fmri_meta_object(meta))
  expect_error(fmrireg:::.validate_fmri_meta_object(list()), "fmri_meta")
  expect_error(
    fmrireg:::.validate_fmri_meta_object(
      structure(list(coefficients = matrix(1, 0, 1)), class = "fmri_meta")
    ),
    "coefficient"
  )
  bad <- meta
  bad$se <- matrix(1, 5, 1)
  expect_error(fmrireg:::.validate_fmri_meta_object(bad), "matching dimensions")

  expect_equal(
    fmrireg:::.fmri_meta_coef_names(meta),
    c("(Intercept)", "groupB")
  )
  anon <- meta
  colnames(anon$coefficients) <- NULL
  expect_equal(fmrireg:::.fmri_meta_coef_names(anon), c("coef1", "coef2"))

  expect_equal(
    fmrireg:::.normalize_fmri_meta_stats(c("beta", "p", "pval")),
    c("beta", "pval")
  )

  idx <- fmrireg:::.select_fmri_meta_coefficients(
    c("(Intercept)", "groupB", "age"),
    coefficients = "groupB",
    coefficient_match = "exact"
  )
  expect_equal(idx, 2L)
  idx_re <- fmrireg:::.select_fmri_meta_coefficients(
    c("(Intercept)", "groupB", "age"),
    coefficients = "group",
    coefficient_match = "regex"
  )
  expect_equal(idx_re, 2L)
  expect_error(
    fmrireg:::.select_fmri_meta_coefficients(
      c("a", "b"), coefficients = "(bad", coefficient_match = "regex"
    ),
    "Invalid coefficient regex"
  )
  expect_error(
    fmrireg:::.select_fmri_meta_coefficients(
      c("a", "b"), coefficients = "zzz", coefficient_match = "exact"
    ),
    "No coefficients matched"
  )

  expect_true(fmrireg:::.fmri_meta_has_heterogeneity(meta))
  expect_false(fmrireg:::.fmri_meta_has_heterogeneity(list()))

  tidy_roi <- tidy(meta, conf.int = TRUE, conf.level = 0.9)
  expect_s3_class(tidy_roi, "tbl_df")
  expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(tidy_roi)))

  tidy_vox <- tidy(toy_meta(n_feat = 6L, roi = FALSE))
  expect_true(all(c("term", "mean_estimate", "n_significant") %in% names(tidy_vox)))
})

test_that("fit_meta_chunk / single / tau2 / huber cover meta internals", {
  set.seed(3)
  n <- 10L
  p <- 4L
  Y <- matrix(rnorm(n * p, mean = 0.3), n, p)
  se_mat <- matrix(runif(n * p, 0.1, 0.3), n, p)
  X <- cbind(1, rnorm(n))

  chunk <- fmrireg:::fit_meta_chunk(
    list(beta = Y, se = se_mat),
    X = X, method = "fe", robust = "none",
    weights = "ivw", weights_custom = NULL, n_threads = 1
  )
  expect_equal(dim(chunk$coefficients), c(p, 2L))
  expect_equal(dim(chunk$se), c(p, 2L))

  chunk_eq <- fmrireg:::fit_meta_chunk(
    list(beta = Y, se = se_mat),
    X = X, method = "dl", robust = "huber",
    weights = "equal", weights_custom = NULL, n_threads = 1
  )
  expect_equal(length(chunk_eq$tau2), p)

  Cmat <- matrix(c(0, 1), nrow = 2)
  chunk_c <- fmrireg:::fit_meta_chunk(
    list(beta = Y, se = se_mat),
    X = X, method = "fe", robust = "none",
    weights = "ivw", weights_custom = NULL, n_threads = 1,
    contrasts = Cmat, return_cov = "tri"
  )
  expect_equal(dim(chunk_c$con_est), c(1L, p))
  expect_true(!is.null(chunk_c$cov_tri))

  # t-only combine path
  tmat <- matrix(rnorm(n * p), n, p)
  attr(tmat, "df") <- 20
  comb <- fmrireg:::fit_meta_chunk(
    list(t = tmat),
    X = matrix(1, n, 1), method = "fe", robust = "none",
    weights = "equal", weights_custom = NULL, n_threads = 1,
    combine = "stouffer"
  )
  expect_equal(nrow(comb$coefficients), p)

  expect_error(
    fmrireg:::fit_meta_chunk(
      list(t = tmat),
      X = X, method = "fe", robust = "none",
      weights = "equal", weights_custom = NULL, n_threads = 1,
      combine = "stouffer"
    ),
    "intercept-only"
  )
  expect_error(
    fmrireg:::fit_meta_chunk(
      list(t = matrix(1, n, p)),
      X = matrix(1, n, 1), method = "fe", robust = "none",
      weights = "equal", weights_custom = NULL, n_threads = 1,
      combine = "stouffer"
    ),
    "without 'df'"
  )
  expect_error(
    fmrireg:::fit_meta_chunk(
      list(t = tmat),
      X = matrix(1, n, 1), method = "fe", robust = "none",
      weights = "equal", weights_custom = NULL, n_threads = 1
    ),
    "combine"
  )
  expect_error(
    fmrireg:::fit_meta_chunk(
      list(),
      X = X, method = "fe", robust = "none",
      weights = "ivw", weights_custom = NULL, n_threads = 1
    ),
    "No usable data"
  )

  single <- fmrireg:::fit_meta_single(
    y = Y[, 1], se = se_mat[, 1], X = X,
    method = "pm", robust = "none", weights = "ivw", weights_custom = NULL
  )
  expect_equal(length(single$coefficients), 2L)
  expect_true(is.finite(single$tau2))

  tau_dl <- fmrireg:::estimate_tau2(Y[, 1], se_mat[, 1], X, method = "dl")
  tau_pm <- fmrireg:::estimate_tau2(Y[, 1], se_mat[, 1], X, method = "pm")
  tau_reml <- fmrireg:::estimate_tau2(Y[, 1], se_mat[, 1], X, method = "reml")
  expect_true(all(c(tau_dl, tau_pm, tau_reml) >= 0))

  w0 <- 1 / se_mat[, 1]^2
  w_h <- fmrireg:::apply_huber_weights(Y[, 1], se_mat[, 1], X, w0)
  expect_equal(length(w_h), n)
  expect_true(all(w_h > 0))

  # parse_meta_formula with synthetic group data
  gd <- structure(
    list(
      subjects = paste0("s", 1:n),
      covariates = data.frame(group = factor(rep(c("A", "B"), length.out = n)))
    ),
    class = c("group_data_csv", "group_data")
  )
  parsed <- fmrireg:::parse_meta_formula(~ group, gd)
  expect_equal(nrow(parsed$X), n)
  expect_true("(Intercept)" %in% colnames(parsed$X))

  gd2 <- structure(list(subjects = paste0("s", 1:5), covariates = NULL),
                   class = c("group_data_csv", "group_data"))
  parsed2 <- fmrireg:::parse_meta_formula(~ 1, gd2)
  expect_equal(parsed2$X, matrix(1, 5, 1, dimnames = list(NULL, "(Intercept)")))
})
