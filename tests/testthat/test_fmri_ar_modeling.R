context("fmri_ar_modeling utilities")

# ensure deterministic results
options(mc.cores=1)

whiten_R <- function(M, phi, exact_first = FALSE) {
  p <- length(phi)
  nT <- nrow(M)
  nC <- ncol(M)
  res <- M
  for (c in seq_len(nC)) {
    prev <- rep(0, p)
    for (t in seq_len(nT)) {
      orig <- M[t, c]
      val <- orig - sum(phi * prev)
      if (t == 1 && exact_first && p == 1) {
        val <- val * sqrt(1 - phi^2)
      }
      if (p > 0) {
        if (p > 1) prev[p:2] <- prev[(p-1):1]
        prev[1] <- orig
      }
      res[t, c] <- val
    }
  }
  res
}

# estimate_ar_parameters should mirror .estimate_ar

test_that("estimate_ar_parameters recovers AR coefficients", {
  set.seed(123)
  phi1 <- 0.6
  x1 <- as.numeric(arima.sim(model = list(ar = phi1), n = 500))
  est1 <- estimate_ar_parameters(x1, 1)
  expect_equal(as.numeric(est1), phi1, tolerance = 0.1)

  phi2 <- c(0.5, -0.25)
  x2 <- as.numeric(arima.sim(model = list(ar = phi2), n = 1000))
  est2 <- estimate_ar_parameters(x2, 2)
  expect_equal(as.numeric(est2), phi2, tolerance = 0.1)
})

# ar_whiten_transform should return whitened matrices

test_that("ar_whiten_transform whitens X and Y", {
  Y <- matrix(1:4, ncol = 1)
  X <- matrix(4:1, ncol = 1)
  phi <- 0.5
  refY <- whiten_R(Y, phi)
  refX <- whiten_R(X, phi)
  res <- ar_whiten_transform(X, Y, phi)
  expect_equal(res$Y, refY)
  expect_equal(res$X, refX)
})

# ar_whiten_transform error on NA

test_that("ar_whiten_transform errors with NA", {
  Y <- matrix(rnorm(4), ncol = 1)
  X <- matrix(rnorm(4), ncol = 1)
  Y[2, 1] <- NA
  expect_error(ar_whiten_transform(X, Y, 0.3), "NA")
})
