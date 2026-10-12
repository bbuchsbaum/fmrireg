# Review 2026-10, pass 3, stream R3-2: exported API surface.
#
# R3-2: `.resample_param` and `.fmri_lm_formula_legacy` were exported
#   dot-functions and `estimate()` an exported stub that always errored.
#   fmrireg defined its own `tidy` generic, so `generics::tidy()` /
#   `broom::tidy()` never reached `tidy.fmri_lm`; it now re-exports
#   `generics::tidy`. One `%||%` definition.

test_that("R3-2: internal dot-functions and the estimate() stub are not exported", {
  exports <- getNamespaceExports("fmrireg")
  expect_false(".resample_param" %in% exports)
  expect_false(".fmri_lm_formula_legacy" %in% exports)
  expect_false("estimate" %in% exports)
  expect_false(exists("estimate", envir = asNamespace("fmrireg"), inherits = FALSE))
  # still available internally for their callers
  expect_true(is.function(fmrireg:::.resample_param))
  expect_true(is.function(fmrireg:::.fmri_lm_formula_legacy))
})

test_that("R3-2: tidy is generics::tidy, so generics/broom dispatch to fmrireg methods", {
  expect_true("tidy" %in% getNamespaceExports("fmrireg"))
  expect_identical(fmrireg::tidy, generics::tidy)

  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  ours <- fmrireg::tidy(fit)
  via_generics <- generics::tidy(fit)
  expect_s3_class(via_generics, "tbl_df")
  expect_identical(via_generics, ours)
  # tidy.fmri_lm's own layout (one row per voxel x term), not a fallback
  expect_true(all(c("voxel", "term", "estimate", "std_error", "statistic",
                    "p_value") %in% names(via_generics)))
  expect_equal(nrow(via_generics), nrow(coef(fit, include_baseline = TRUE)) *
                 ncol(coef(fit, include_baseline = TRUE)))

  skip_if_not_installed("broom")
  expect_identical(broom::tidy(fit), ours)
})

test_that("R3-2: a single %||% definition, in utils_operators.R", {
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package sources not available")
  defs <- unlist(lapply(list.files(r_dir, "\\.R$", full.names = TRUE), function(f) {
    hits <- grep("^`%\\|\\|%`\\s*<-", readLines(f, warn = FALSE))
    rep(basename(f), length(hits))
  }))
  expect_identical(defs, "utils_operators.R")
  `%||%` <- fmrireg:::`%||%`
  expect_identical(NULL %||% 2, 2)
  expect_identical(1 %||% 2, 1)
})
