test_that("fmri_meta.gds matches legacy CSV path on demo data", {
  skip_if_not_installed("fmrigds")

  # Legacy CSV group_data path
  gd_csv <- fmrireg:::.demo_group_data_csv()
  fit_legacy <- fmri_meta(gd_csv, formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)

  # Build equivalent gds via group_data wrapper
  df <- gd_csv$data
  gd_gds <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )

  fit_gds <- fmri_meta(gd_gds, formula = ~ 1, method = "fe", robust = "none", verbose = FALSE)

  expect_s3_class(fit_gds, "fmri_meta")
  expect_equal(dim(fit_gds$coefficients), dim(fit_legacy$coefficients))
  expect_equal(dim(fit_gds$se), dim(fit_legacy$se))

  # Numeric parity within tolerance
  expect_equal(fit_gds$coefficients, fit_legacy$coefficients, tolerance = 1e-8)
  expect_equal(fit_gds$se,           fit_legacy$se,           tolerance = 1e-6)
})

test_that("fmri_meta.gds returns covariance triangles when requested", {
  skip_if_not_installed("fmrigds")

  gd_csv <- fmrireg:::.demo_group_data_csv()
  df <- gd_csv$data
  gd_gds <- group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )

  fit_cov <- fmri_meta(gd_gds, formula = ~ 1, method = "pm", return_cov = "tri", verbose = FALSE)
  expect_true(!is.null(fit_cov$cov))
  expect_equal(fit_cov$cov$type, "tri")
  # Diagonal of Var(beta) equals se^2 for intercept-only model
  se2 <- fit_cov$se^2
  # cov$tri is tsize x P; for intercept-only, tsize = 1 and equals Var(beta)
  expect_equal(as.numeric(fit_cov$cov$tri), as.numeric(se2), tolerance = 1e-6)
})
