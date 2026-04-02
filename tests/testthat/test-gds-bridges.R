test_that("group_data_gds bridges work", {
  skip_if_not_installed("fmrigds")
  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )
  expect_equal(n_subjects(gd), length(unique(df$subject)))
  expect_equal(get_subjects(gd), unique(df$subject))
  expect_true(is.data.frame(get_covariates(gd)))
})

test_that("group_data_gds preserves sample metadata for fmrireg helpers", {
  skip_if_not_installed("fmrigds")
  df <- data.frame(
    subject = rep(paste0("s", 1:4), each = 3),
    roi = rep(c("ROI1", "ROI2", "ROI3"), times = 4),
    contrast = "A_vs_B",
    beta = rnorm(12),
    se = rep(0.2, 12)
  )

  gd <- group_data(
    df,
    format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    feature_group = c("g1", "g1", "g2")
  )

  expect_equal(attr(gd, "fmrireg_sample_labels", exact = TRUE), c("ROI1", "ROI2", "ROI3"))
  expect_equal(as.character(attr(gd, "fmrireg_feature_group", exact = TRUE)), c("g1", "g1", "g2"))
  expect_equal(fmrireg:::.fmri_ttest_sample_labels(gd, NULL), c("ROI1", "ROI2", "ROI3"))
  expect_equal(as.character(fmrireg:::.fmri_ttest_feature_group(gd, NULL)), c("g1", "g1", "g2"))
})
