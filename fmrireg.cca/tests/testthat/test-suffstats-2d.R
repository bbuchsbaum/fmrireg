test_that("2D suffstats equals explicit S cross-products", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)
  set.seed(1)
  nx <- 12; ny <- 10; T <- 40
  vols <- replicate(T, matrix(rnorm(nx*ny), nx, ny), simplify = FALSE)
  spacing <- c(2,2); fwhm <- 6
  V <- nx*ny; mask_lin <- seq_len(V); mask_z <- rep(1L, V)
  w_dir   <- matrix(c(0.5, 0.3, 0.2), V, 3, byrow = TRUE)
  w_step2 <- matrix(c(0.7, 0.3), V, 2, byrow = TRUE)
  X <- matrix(rnorm(T), T, 1)
  acc <- fmrireg.cca:::friman_pass2_xts_sts_2d(vols, spacing, fwhm, X, w_dir, w_step2, mask_lin, mask_z)
  S <- fmrireg.cca:::friman_apply_series_2d(vols, spacing, fwhm, w_dir, w_step2, mask_lin, mask_z)
  XtX <- crossprod(X); XtS <- crossprod(X, S); StS <- colSums(S^2)
  expect_equal(acc$XtS, XtS, tolerance = 1e-8)
  expect_equal(acc$StS, StS, tolerance = 1e-8)
})

