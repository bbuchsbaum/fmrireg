context("AR coefficient estimation")

options(mc.cores=1)

library(testthat)

# .estimate_ar should recover AR coefficients from simulated data

test_that(".estimate_ar recovers AR(1) coefficient", {
  set.seed(123)
  phi <- 0.6
  x <- as.numeric(arima.sim(model = list(ar = phi), n = 500))
  est <- .estimate_ar(x, 1)
  expect_equal(as.numeric(est), phi, tolerance = 0.05)
})

test_that(".estimate_ar recovers AR(2) coefficients", {
  set.seed(123)
  phi <- c(0.5, -0.25)
  x <- as.numeric(arima.sim(model = list(ar = phi), n = 1000))
  est <- .estimate_ar(x, 2)
  expect_equal(as.numeric(est), phi, tolerance = 0.05)
})
