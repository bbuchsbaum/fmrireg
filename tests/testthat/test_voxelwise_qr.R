context("Voxelwise QR contrast engine")

library(fmrireg)
set.seed(123)

test_that("voxelwise QR and XtXinv methods give same results", {
  n <- 20
  p <- 3
  V <- 5

  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- matrix(rnorm(n * V), n, V)

  proj <- fmrireg:::.fast_preproject(X)
  fit <- fmrireg:::.fast_lm_matrix(X, Y, proj)

  Betas <- fit$betas
  sigma_vec <- sqrt(fit$sigma2)
  qr_obj <- proj$qr

  XtXinv_list <- replicate(V, proj$XtXinv, simplify = FALSE)
  qr_list <- replicate(V, qr_obj, simplify = FALSE)

  # Create contrast properly - colind should indicate which columns the contrast applies to
  conlist <- list(A = 1)  # Simple contrast for coefficient 2
  attr(conlist$A, "colind") <- 2L  # This contrast applies to column 2
  fconlist <- list()

  res_xtx <- fmrireg:::fit_lm_contrasts_voxelwise(Betas, sigma_vec^2, XtXinv_list,
                                         conlist, fconlist, proj$dfres)
  res_qr <- fmrireg:::fit_lm_contrasts_voxelwise_qr(Betas, qr_list, sigma_vec,
                                          conlist, fconlist, proj$dfres)

  expect_equal(res_qr$A$data[[1]], res_xtx$A$data[[1]], tolerance = 1e-2)
})
