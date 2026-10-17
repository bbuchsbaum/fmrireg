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
  # Walk the parsed calls rather than grepping deparsed text: coverage
  # instrumentation wraps each expression in covr:::count(), which is not
  # code of .onLoad's own.
  private_refs <- function(x) {
    if (!is.call(x)) return(character(0))
    fn <- x[[1]]
    if (identical(fn, as.name(":::"))) {
      pkg <- as.character(x[[2]])
      return(if (identical(pkg, "covr")) character(0) else paste0(pkg, ":::", as.character(x[[3]])))
    }
    hit <- if (is.name(fn) && as.character(fn) %in% c("getFromNamespace", "asNamespace")) {
      as.character(fn)
    } else if (is.call(fn) && identical(fn[[1]], as.name("::")) &&
               as.character(fn[[3]]) %in% c("getFromNamespace", "asNamespace")) {
      as.character(fn[[3]])
    }
    c(hit, unlist(lapply(as.list(x), private_refs), use.names = FALSE))
  }
  expect_identical(private_refs(body(fmrireg:::.onLoad)), character(0))

  # the walker does catch the forbidden forms
  probe <- quote({ a <- fmridesign:::foo; b <- utils::getFromNamespace("x", "y"); asNamespace("z") })
  expect_identical(private_refs(probe),
                   c("fmridesign:::foo", "getFromNamespace", "asNamespace"))
  expect_identical(private_refs(quote(covr:::count("id", f(x)))), character(0))
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
