# tests/testthat/test-ar-whitening.R

test_that("apply_ar_filter collapses AR(1) & AR(2) autocorrelation", {
  skip_on_cran()

  # Helpers --------------------------------------------------------------
  sim_ar1 <- function(T, K, rho, sd = 1) {
    M <- matrix(0, T, K)
    M[1, ] <- rnorm(K, sd = sd / sqrt(1 - rho^2))
    for (t in 2:T) M[t, ] <- rho * M[t - 1, ] + rnorm(K, sd = sd)
    M
  }
  sim_ar2 <- function(T, K, phi1, phi2, sd = 1) {
    M <- matrix(0, T, K)
    M[1, ] <- rnorm(K, sd = sd)
    if (T > 1) M[2, ] <- phi1 * M[1, ] + rnorm(K, sd = sd)
    if (T > 2) for (t in 3:T) M[t, ] <- phi1 * M[t - 1, ] + phi2 * M[t - 2, ] + rnorm(K, sd = sd)
    M
  }
  lagk <- function(x, k = 1L) {
    x <- x - mean(x)
    num <- sum(head(x, -k) * tail(x, -k))
    den <- sum(x * x)
    if (den == 0) return(0)
    num / den
  }

  # Parameters -----------------------------------------------------------
  set.seed(123)
  Tlen <- 400L
  K    <- 64L

  # --- AR(1) case ------------------------------------------------------
  rho <- 0.7
  M1  <- sim_ar1(Tlen, K, rho = rho, sd = 1.0)
  a1_before <- mean(apply(M1, 2, lagk, k = 1L))

  # Create dummy design matrix
  X1 <- matrix(1, nrow = Tlen, ncol = 1)
  
  # Apply AR whitening
  result1 <- fmrireg:::ar_whiten_inplace(M1, X1, phi_coeffs = rho, exact_first_ar1 = FALSE)
  M1w <- result1$Y
  a1_after  <- mean(apply(M1w, 2, lagk, k = 1L))

  # Diagnostics: near-zero and much smaller than before
  expect_lt(abs(a1_after), 0.08)                # ~1/sqrt(T) margin with slack
  expect_lt(abs(a1_after), abs(a1_before) / 5)  # big shrink vs. raw

  # --- AR(2) case ------------------------------------------------------
  phi1 <- 0.5; phi2 <- -0.2
  M2   <- sim_ar2(Tlen, K, phi1 = phi1, phi2 = phi2, sd = 1.0)
  a1b  <- mean(apply(M2, 2, lagk, k = 1L))
  a2b  <- mean(apply(M2, 2, lagk, k = 2L))

  # Create dummy design matrix
  X2 <- matrix(1, nrow = Tlen, ncol = 1)
  
  # Apply AR whitening
  result2 <- fmrireg:::ar_whiten_inplace(M2, X2, phi_coeffs = c(phi1, phi2), exact_first_ar1 = FALSE)
  M2w <- result2$Y
  a1a  <- mean(apply(M2w, 2, lagk, k = 1L))
  a2a  <- mean(apply(M2w, 2, lagk, k = 2L))

  # Diagnostics: both lags shrink a lot and are small in absolute value
  expect_lt(abs(a1a), max(0.12, abs(a1b) / 3))
  expect_lt(abs(a2a), max(0.12, abs(a2b) / 3))
})