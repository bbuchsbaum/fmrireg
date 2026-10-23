# Simplification 2026-10, legacy group-data deletion prep: behaviour the
# legacy csv reader provided that gds group data now provide too, so deleting
# the legacy stack loses nothing. Each test failed before the prep.

skip_if_not_installed("fmrigds")

prep_table <- function() {
  data.frame(
    subject = rep(sprintf("s%02d", 1:6), each = 2),
    roi = rep(c("ACC", "PFC"), times = 6),
    beta = c(0.3, 0.1, 0.5, 0.2, 0.2, 0.4, 0.6, 0.1, 0.4, 0.3, 0.1, 0.2),
    se = rep(c(0.1, 0.2), times = 6)
  )
}

test_that("fmri_meta on ROI-table gds data is an ROI fit (tidy per ROI, prints ROIs)", {
  gd <- group_data(prep_table(), format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject", roi_col = "roi")
  fit <- fmri_meta(gd, ~ 1, method = "fe", verbose = FALSE)
  expect_s3_class(fit, "fmri_meta_roi")
  expect_equal(fit$n_rois, 2)
  td <- tidy(fit)
  expect_identical(td$roi, c("ACC", "PFC"))
  for (r in c("ACC", "PFC")) {
    d <- prep_table()[prep_table()$roi == r, ]
    w <- 1 / d$se^2
    expect_equal(td$estimate[td$roi == r], sum(w * d$beta) / sum(w), tolerance = 1e-12)
    expect_equal(td$std.error[td$roi == r], 1 / sqrt(sum(w)), tolerance = 1e-12)
  }
  expect_output(print(fit), "ROIs analyzed: 2")
})

test_that("a csv table without an ROI column is one ROI, 'overall'", {
  d <- prep_table()[prep_table()$roi == "ACC", c("subject", "beta", "se")]
  gd <- group_data(d, format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject")
  fit <- fmri_meta(gd, ~ 1, method = "fe", verbose = FALSE)
  expect_identical(fit$roi_names, "overall")
  w <- 1 / d$se^2
  expect_equal(unname(fit$coefficients[1, 1]), sum(w * d$beta) / sum(w), tolerance = 1e-12)
})

test_that("summary() on gds group data uses the gds method, not the legacy one", {
  gd <- group_data(prep_table(), format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject", roi_col = "roi")
  out <- capture.output(summary(gd))
  expect_true(any(grepl("group_data (backed by gds)", out, fixed = TRUE)))
  expect_false(any(grepl("^Format:", out)))
})
