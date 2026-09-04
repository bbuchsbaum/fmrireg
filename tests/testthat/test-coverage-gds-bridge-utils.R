# Coverage for gds_bridge_utils.R and group_data_gds annotation/print helpers.

test_that(".gds_coerce_matrix handles matrix/df/array/vector/null", {
  expect_null(fmrireg:::.gds_coerce_matrix(NULL))

  m <- matrix(1:6, 2, 3)
  expect_identical(fmrireg:::.gds_coerce_matrix(m), m)

  df <- data.frame(a = 1:3, b = 4:6)
  expect_equal(fmrireg:::.gds_coerce_matrix(df), as.matrix(df))

  v <- 1:4
  expect_equal(fmrireg:::.gds_coerce_matrix(v), matrix(1:4, ncol = 1))

  a1 <- array(1:5, dim = 5)
  expect_equal(dim(fmrireg:::.gds_coerce_matrix(a1)), c(5L, 1L))

  a2 <- array(1:12, dim = c(3, 2, 2))
  coerced <- fmrireg:::.gds_coerce_matrix(a2)
  expect_equal(dim(coerced), c(3L, 4L))

  # Degenerate array with NULL dim via as.vector path: structure without dim
  bare <- structure(1:3, class = "array")
  # If dim is NULL after as.array quirks, coerce still returns something usable
  out_bare <- fmrireg:::.gds_coerce_matrix(bare)
  expect_true(is.matrix(out_bare) || is.atomic(out_bare))
})

test_that(".gds_terms/features helpers and safe assay fallbacks", {
  expect_null(fmrireg:::.gds_terms_by_features(NULL))
  expect_null(fmrireg:::.gds_features_by_terms(NULL))

  m <- matrix(1:6, 2, 3, dimnames = list(c("t1", "t2"), c("f1", "f2", "f3")))
  expect_equal(fmrireg:::.gds_terms_by_features(m), t(m))
  expect_equal(fmrireg:::.gds_features_by_terms(m), m)

  # Safe assay helpers return NULL/character(0) without crashing on non-gds
  expect_null(fmrireg:::.gds_safe_assay(NULL, "beta"))
  expect_equal(fmrireg:::.gds_safe_assay_names(NULL), character(0))
  expect_null(fmrireg:::.gds_safe_assay(list(), "beta"))
  expect_equal(fmrireg:::.gds_safe_assay_names(list()), character(0))
})

test_that(".gds_safe_model_matrix falls back to intercept design", {
  # Non-gds object: tryCatch fails, intercept formula yields fallback
  mm <- fmrireg:::.gds_safe_model_matrix(list(), ~ 1, fallback_n = 5L)
  expect_equal(dim(mm), c(5L, 1L))
  expect_equal(colnames(mm), "(Intercept)")

  expect_null(fmrireg:::.gds_safe_model_matrix(list(), ~ 1, fallback_n = NULL))
  expect_null(fmrireg:::.gds_safe_model_matrix(list(), ~ age, fallback_n = 5L))
})

test_that(".annotate_group_data_gds and print/accessors work on stubs", {
  stub <- structure(list(ok = TRUE), class = "dummy_gds")
  annotated <- fmrireg:::.annotate_group_data_gds(
    stub,
    sample_labels = c("roi1", "roi2"),
    feature_group = c("g1", "g2")
  )
  expect_true(inherits(annotated, "group_data_gds"))
  expect_true(inherits(annotated, "group_data"))
  expect_equal(attr(annotated, "fmrireg_sample_labels"), c("roi1", "roi2"))
  expect_equal(attr(annotated, "fmrireg_feature_group"), c("g1", "g2"))

  # Print should not error even without fmrigds assays
  expect_output(
    fmrireg:::.print_group_data_gds(annotated),
    "group_data"
  )

  # Feature/sample label helpers prefer attributes, then block metadata
  expect_equal(
    fmrireg:::.fmri_ttest_feature_group(annotated, NULL),
    c("g1", "g2")
  )
  expect_equal(
    fmrireg:::.fmri_ttest_sample_labels(annotated, NULL),
    c("roi1", "roi2")
  )
  expect_equal(
    fmrireg:::.fmri_ttest_feature_group(
      list(),
      block = list(feature = list(group = c("a", "b")))
    ),
    c("a", "b")
  )
  expect_equal(
    fmrireg:::.fmri_ttest_sample_labels(
      list(),
      block = list(feature = list(label = c("x", "y")))
    ),
    c("x", "y")
  )
  expect_null(fmrireg:::.fmri_ttest_feature_group(list(), NULL))
  expect_null(fmrireg:::.fmri_ttest_sample_labels(list(), NULL))
})

test_that(".group_data_gds builds from tabular data when fmrigds is present", {
  skip_if_not_installed("fmrigds")
  df <- data.frame(
    subject = rep(paste0("s", 1:3), each = 2),
    roi = rep(c("ROI1", "ROI2"), times = 3),
    contrast = "A_vs_B",
    age = rep(c(20, 30, 40), each = 2),
    beta = c(0.1, 0.2, 0.15, 0.25, 0.2, 0.3),
    se = rep(0.1, 6)
  )

  gd <- fmrireg:::.group_data_gds(
    df,
    format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age",
    feature_group = c("g1", "g2")
  )
  expect_s3_class(gd, "group_data_gds")
  expect_equal(attr(gd, "fmrireg_sample_labels"), c("ROI1", "ROI2"))
  expect_equal(as.character(attr(gd, "fmrireg_feature_group")), c("g1", "g2"))

  expect_error(
    fmrireg:::.group_data_gds(
      df,
      format = "csv",
      effect_cols = c(beta = "beta", se = "se"),
      subject_col = "subject",
      roi_col = "roi",
      contrast_col = "contrast",
      feature_group = c("only_one")
    ),
    "feature_group must have length"
  )

  # Accessors
  ns <- fmrireg:::.n_subjects_group_data_gds(gd)
  expect_equal(ns, 3L)
  subs <- fmrireg:::.get_subjects_group_data_gds(gd)
  expect_equal(length(subs), 3L)
  covs <- fmrireg:::.get_covariates_group_data_gds(gd)
  expect_true(is.data.frame(covs) || is.null(covs) || is.matrix(covs))

  expect_output(fmrireg:::.print_group_data_gds(gd), "group_data|assays|subjects")
})
