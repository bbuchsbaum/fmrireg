# Review 2026-10, pass 3, stream R3-5: bootstrap_glm_inference() was internal,
# never called by fmri_lm(), exercised only by tests, and statistically
# defective (the block bootstrap concatenated blocks across runs, the wild
# bootstrap destroyed autocorrelation, the BCa interval returned NaN). It is
# removed with its helpers, as is compute_contrast(), whose only callers were
# the bootstrap.

test_that("R3-5: the unwired bootstrap and its single-use contrast helper are gone", {
  ns <- asNamespace("fmrireg")
  removed <- c("bootstrap_glm_inference", "create_bootstrap_blocks",
               "bootstrap_residual", "bootstrap_case", "bootstrap_wild",
               "bootstrap_hypothesis_test", "compute_bca_ci", "compute_contrast")
  for (f in removed) {
    expect_false(exists(f, envir = ns, inherits = FALSE), info = f)
  }
  # the progress package was only imported for the bootstrap's progress bar
  expect_false("progress" %in% names(getNamespaceImports("fmrireg")))
})
