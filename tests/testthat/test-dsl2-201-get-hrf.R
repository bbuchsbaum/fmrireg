library(testthat)

GH <- fmrireg:::get_hrf_from_dsl

make_cfg <- function(hrf_list) {
  list(hrfs = hrf_list)
}

test_that("SPMCanonical maps to HRF_SPMG1", {
  cfg <- make_cfg(list(canonical = list(type = "SPMCanonical")))
  h <- GH("canonical", cfg)
  expect_s3_class(h, "HRF")
  expect_equal(nbasis(h), 1L)
})

test_that("SPMCanonicalDerivs derivatives handled", {
  cfg <- make_cfg(list(cand = list(type="SPMCanonicalDerivs", derivatives=list("Temporal"))))
  h <- GH("cand", cfg)
  expect_equal(nbasis(h), 2L)

  cfg2 <- make_cfg(list(cand = list(type="SPMCanonicalDerivs", derivatives=list("Temporal","Dispersion"))))
  h2 <- GH("cand", cfg2)
  expect_equal(nbasis(h2), 3L)
})

test_that("GammaFunction parameters applied", {
  cfg <- make_cfg(list(gam = list(type="GammaFunction", parameters=list(shape=7, rate=0.5))))
  h <- GH("gam", cfg)
  p <- attr(h, "params")
  expect_equal(p$shape, 7)
  expect_equal(p$rate, 0.5)
})

test_that("CustomR loads function", {
  myfun <- function(t) t*0
  cfg <- make_cfg(list(custom=list(type="CustomR", definition="myfun", lag=2)))
  h <- GH("custom", cfg)
  expect_true(inherits(h, "HRF"))
  expect_equal(attr(h, "params")$.lag, 2)
})

test_that("CustomR evaluates code string", {
  cfg <- make_cfg(list(custom=list(type="CustomR", definition="function(t) t*0")))
  h <- GH("custom", cfg)
  expect_s3_class(h, "HRF")
  expect_equal(as.numeric(h(0)), 0)
})

