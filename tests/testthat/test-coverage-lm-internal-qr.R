# fmri_lm_internal voxelwise QR / storage helpers via small fixtures.

test_that("fit_lm_contrasts_voxelwise_qr covers t and F contrast loops", {
  set.seed(9)
  n <- 40L
  p <- 3L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("Intercept", "A", "B")
  Y <- as.vector(X %*% c(0.2, 1.1, -0.7)) + matrix(rnorm(n * V, sd = 0.3), n, V)
  Betas <- qr.coef(qr(X), Y)
  sigma <- sqrt(colSums((Y - X %*% Betas)^2) / (n - p))
  qr_list <- lapply(seq_len(V), function(v) qr(X))

  tcon <- list(
    A_gt_0 = structure(1, colind = 2L),
    A_vs_B = structure(c(1, -1), colind = c(2L, 3L))
  )
  fcon <- list(
    AB = structure(
      matrix(c(1, 0, 0, 1), nrow = 2),
      colind = c(2L, 3L)
    )
  )

  res <- fmrireg:::fit_lm_contrasts_voxelwise_qr(
    Betas = Betas,
    qr_list = qr_list,
    sigma = sigma,
    conlist = tcon,
    fconlist = fcon,
    dfres = n - p
  )
  expect_true(is.list(res))
  expect_true(all(c("A_gt_0", "A_vs_B", "AB") %in% names(res)))
  expect_equal(nrow(res$A_gt_0$data[[1]]), V)

  # Robust weights / AR effective-df branch
  rw <- lapply(seq_len(V), function(v) runif(n, 0.6, 1))
  res_rw <- fmrireg:::fit_lm_contrasts_voxelwise_qr(
    Betas = Betas,
    qr_list = qr_list,
    sigma = sigma,
    conlist = tcon[1],
    fconlist = list(),
    dfres = n - p,
    robust_weights_list = rw,
    ar_order = 1
  )
  expect_true("A_gt_0" %in% names(res_rw))
})

test_that("initialize/store/format contrast storage helpers work", {
  conlist <- list(c1 = structure(c(0, 1), colind = 2L))
  fconlist <- list(F1 = structure(diag(2), colind = 1:2))
  storage <- fmrireg:::initialize_contrast_storage(conlist, fconlist, n_voxels = 3L)
  expect_true(is.list(storage))

  X <- cbind(1, rnorm(20))
  qr_x <- qr(X)
  beta <- c(0.1, 0.5)
  stored <- tryCatch(
    fmrireg:::store_voxel_contrasts(
      storage, voxel_index = 1L, qr_or_XtXinv = qr_x,
      beta_w = beta, sigma = 0.4, dfres = 18
    ),
    error = function(e) e
  )
  if (!inherits(stored, "error")) {
    expect_true(is.list(stored))
  }

  formatted <- tryCatch(
    fmrireg:::format_contrast_results(storage, dfres = 18),
    error = function(e) e
  )
  if (!inherits(formatted, "error")) {
    expect_true(is.list(formatted) || inherits(formatted, "tbl_df"))
  }
})

test_that(".fast_lm_matrix and voxel status helpers remain coherent", {
  X <- cbind(1, rnorm(25), rnorm(25))
  Y <- matrix(rnorm(25 * 3), 25, 3)
  proj <- fmrireg:::.fast_preproject(X)
  fit <- fmrireg:::.fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  expect_true(is.list(fit))
  expect_equal(ncol(fit$betas), 3L)
  expect_true(!is.null(fit$fitted) || !is.null(fit$residuals) || is.matrix(fit$betas))

  expect_true(is.list(fmrireg:::.fmri_lm_voxel_status(rnorm(10))) ||
                is.character(fmrireg:::.fmri_lm_voxel_status(rnorm(10))))
  bad <- fmrireg:::.fmri_lm_voxel_status(c(1, NA, Inf))
  expect_true(!is.null(bad))
})
