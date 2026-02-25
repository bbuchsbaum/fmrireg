context("C++ AR whitening")

whiten_R <- function(M, phi, exact_first_ar1 = FALSE) {
  p <- length(phi)
  nT <- nrow(M)
  nC <- ncol(M)
  res <- M
  for (c in seq_len(nC)) {
    prev <- rep(0, p)
    for (t in seq_len(nT)) {
      orig <- M[t, c]
      val <- orig - sum(phi * prev)
      if (t == 1 && exact_first_ar1 && p == 1) {
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

test_that("ar_whiten_inplace AR1", {
  Y <- matrix(1:4, ncol = 1)
  X <- matrix(4:1, ncol = 1)
  phi <- 0.5
  Yref <- whiten_R(Y, phi)
  Xref <- whiten_R(X, phi)
  result <- ar_whiten_inplace(Y, X, phi)
  expect_equal(result$Y, Yref)
  expect_equal(result$X, Xref)
})

test_that("ar_whiten_inplace AR2", {
  Y <- matrix(1:6, ncol = 1)
  X <- matrix(6:1, ncol = 1)
  phi <- c(0.6, -0.3)
  Yref <- whiten_R(Y, phi)
  Xref <- whiten_R(X, phi)
  result <- ar_whiten_inplace(Y, X, phi)
  expect_equal(result$Y, Yref)
  expect_equal(result$X, Xref)
})

test_that("ar_whiten_inplace exact first AR1", {
  Y <- matrix(rnorm(5), ncol = 1)
  X <- matrix(rnorm(5), ncol = 1)
  phi <- 0.4
  Yref <- whiten_R(Y, phi, TRUE)
  Xref <- whiten_R(X, phi, TRUE)
  result <- ar_whiten_inplace(Y, X, phi, TRUE)
  expect_equal(result$Y, Yref)
  expect_equal(result$X, Xref)
})

test_that("ar_whiten_transform errors with NA input", {
  # Test the R wrapper which does check for NA
  Y <- matrix(rnorm(4), ncol = 1)
  X <- matrix(rnorm(4), ncol = 1)
  Y[2, 1] <- NA
  expect_error(ar_whiten_transform(X, Y, 0.3), "NA")

  Y[2, 1] <- 0
  X[3, 1] <- NA
  expect_error(ar_whiten_transform(X, Y, 0.3), "NA")
})
