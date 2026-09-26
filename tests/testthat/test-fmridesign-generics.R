# fmrireg uses fmridesign's generics by plain dispatch: it defines none of
# its own copies and reaches no private fmridesign namespace object.

test_that("fmrireg's design-name generics are fmridesign's", {
  for (g in c("longnames", "shortnames", "construct", "correlation_map",
              "Fcontrasts", "columns", "cells", "conditions")) {
    expect_identical(getExportedValue("fmrireg", g),
                     getExportedValue("fmridesign", g), info = g)
  }
})

test_that("cross-dispatch: fmridesign:: generics reach fmrireg methods", {
  fm <- fmrireg:::.demo_fmri_model()
  expect_s3_class(fmridesign::correlation_map(fm), "ggplot")
  expect_true(is.data.frame(fmridesign::cells(fm)))
  expect_type(fmridesign::conditions(fm), "character")
})

test_that("cross-dispatch: fmrireg:: re-exports reach fmridesign methods", {
  em <- fmrireg:::.demo_event_model()
  expect_equal(fmrireg::longnames(em), fmridesign::conditions(em, style = "canonical"))
  expect_equal(fmrireg::columns(em), colnames(fmridesign::design_matrix(em)))
  expect_s3_class(fmrireg::correlation_map(em), "ggplot")
  expect_type(fmrireg::Fcontrasts(em), "list")
})

test_that(".onLoad touches no private namespace", {
  src <- paste(deparse(body(fmrireg:::.onLoad)), collapse = "\n")
  expect_false(grepl("getFromNamespace|asNamespace|:::", src))
})

test_that("fmrireg registers no S3 methods for fmridesign-owned classes", {
  fd_classes <- c("event_term", "event_seq", "event", "event_model",
                  "convolved_term", "covariate_convolved_term",
                  "covariate_term", "feature_term", "hrfspec",
                  "baseline_model", "baseline_term", "baselinespec",
                  "blockspec", "nuisancespec", "covariatespec", "featurespec",
                  "ParametricBasis", "Poly", "BSpline", "Ident", "Scale",
                  "ScaleWithin", "Standardized", "RobustScale")
  info <- getNamespaceInfo(asNamespace("fmrireg"), "S3methods")
  offending <- info[info[, 2] %in% fd_classes, , drop = FALSE]
  expect_equal(nrow(offending), 0L,
               info = paste(offending[, 1], offending[, 2], sep = ".", collapse = ", "))
})
