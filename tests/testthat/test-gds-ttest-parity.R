test_that("fmri_ttest.gds (meta engine) matches legacy CSV path", {
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

  fit_old <- fmri_ttest(gd_csv, formula = ~ 1, engine = "meta")
  fit_new <- fmri_ttest(gd_gds, formula = ~ 1, engine = "meta")

  expect_equal(dim(fit_new$beta), dim(fit_old$beta))
  expect_equal(dim(fit_new$se),   dim(fit_old$se))

  expect_equal(fit_new$beta, fit_old$beta, tolerance = 1e-8)
  expect_equal(fit_new$se,   fit_old$se,   tolerance = 1e-6)
  expect_equal(fit_new$z,    fit_old$z,    tolerance = 1e-6)
  expect_equal(fit_new$p,    fit_old$p,    tolerance = 1e-6)
})

test_that("fmri_ttest.gds BH/BY match p.adjust", {
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
  fit_new <- fmri_ttest(gd_gds, formula = ~ 1, engine = "meta", mc = "bh")
  # Compare BH per coefficient row
  for (i in seq_len(nrow(fit_new$p))) {
    expect_equal(as.numeric(fit_new$q[i,]), stats::p.adjust(as.numeric(fit_new$p[i,]), method = "BH"))
  }
  fit_new_by <- fmri_ttest(gd_gds, formula = ~ 1, engine = "meta", mc = "by")
  for (i in seq_len(nrow(fit_new_by$p))) {
    expect_equal(as.numeric(fit_new_by$q[i,]), stats::p.adjust(as.numeric(fit_new_by$p[i,]), method = "BY"))
  }
})
