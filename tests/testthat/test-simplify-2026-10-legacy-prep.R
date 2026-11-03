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

# LEGACY-GROUP-DELETE (next minor): checks the deprecation warnings of the legacy entry points.
test_that("every public legacy group-data entry point warns that it is deprecated", {
  gd_csv <- .demo_group_data_csv()  # built with warnings silenced
  withr::local_options(fmrireg.suppress_deprecation = FALSE)
  expect_warning(extract_csv_data(gd_csv), "extract_csv_data\\(\\) is deprecated")
  expect_warning(get_rois(gd_csv), "get_rois\\(\\) is deprecated")
  expect_warning(get_contrasts(gd_csv), "get_contrasts\\(\\) is deprecated")
  expect_warning(
    tryCatch(fmri_ttest(gd_csv, ~ 1, engine = "meta"), error = function(e) NULL),
    "fmri_ttest\\(\\) on a group_data_csv object is deprecated"
  )
  fake_nifti <- structure(list(format = "nifti", subjects = "s1"),
                          class = c("group_data_nifti", "group_data"))
  expect_warning(
    tryCatch(read_nifti_full(fake_nifti), error = function(e) NULL),
    "read_nifti_full\\(\\) is deprecated"
  )
  expect_warning(
    tryCatch(fmri_meta.group_data_nifti(fake_nifti, verbose = FALSE), error = function(e) NULL),
    "group_data_nifti object is deprecated"
  )
  fake_h5 <- structure(list(format = "h5", subjects = "s1", paths = "none.h5"),
                       class = c("group_data_h5", "group_data"))
  expect_warning(
    tryCatch(read_h5_full(fake_h5), error = function(e) NULL),
    "read_h5_full\\(\\) is deprecated"
  )
  # Internal legacy reads stay quiet: one warning per public call, not several.
  expect_silent(.legacy_group_quiet(extract_csv_data(gd_csv)))
})
