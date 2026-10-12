# fmrigds is a Suggests package. group_data() needs it only for the routes
# that build a gds object; the legacy NIfTI t/df and beta/var routes must work
# without it, and the gds routes must fail with a clear message.

test_that("legacy NIfTI group_data routes work when fmrigds is unavailable", {
  skip_if_not_installed("neuroim2")
  testthat::local_mocked_bindings(
    .pkg_available = function(pkg) !identical(pkg, "fmrigds"),
    .package = "fmrireg")
  set.seed(9)
  dims <- c(2, 2, 1)
  sp <- neuroim2::NeuroSpace(dims)
  write_one <- function(vals) {
    p <- tempfile(fileext = ".nii")
    neuroim2::write_vol(neuroim2::NeuroVol(array(vals, dim = dims), sp), p)
    p
  }
  tstat <- c(write_one(rnorm(4)), write_one(rnorm(4)))
  beta <- c(write_one(rnorm(4)), write_one(rnorm(4)))
  vars <- c(write_one(abs(rnorm(4)) + 0.05), write_one(abs(rnorm(4)) + 0.05))
  on.exit(unlink(c(tstat, beta, vars)), add = TRUE)

  gd_t <- group_data(list(t = tstat, df = 20), format = "nifti",
                     subjects = c("s1", "s2"))
  expect_s3_class(gd_t, "group_data")
  gd_v <- group_data(list(beta = beta, var = vars), format = "nifti",
                     subjects = c("s1", "s2"))
  expect_s3_class(gd_v, "group_data")
})

test_that("gds-backed group_data routes fail clearly when fmrigds is unavailable", {
  testthat::local_mocked_bindings(
    .pkg_available = function(pkg) !identical(pkg, "fmrigds"),
    .package = "fmrireg")
  df <- data.frame(subject = sprintf("s%02d", 1:4), roi = "A",
                   beta = rnorm(4), se = 0.5)
  expect_error(
    group_data(df, format = "csv", effect_cols = c(beta = "beta", se = "se"),
               subject_col = "subject", roi_col = "roi"),
    "requires the 'fmrigds' package")
})
