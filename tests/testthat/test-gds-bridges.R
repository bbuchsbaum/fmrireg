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
