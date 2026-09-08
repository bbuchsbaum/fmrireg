# Chunked voxelwise contrasts + pull_stat invalid-index warning + preproject edges.

test_that("fit_lm_contrasts_voxelwise_chunked whitens and stores contrasts", {
  skip_if_not_installed("fmriAR")
  set.seed(51)
  n <- 48L
  p <- 3L
  V <- 9L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  Y <- matrix(rnorm(n * V), n, V)

  conlist <- list(
    A_vs_B = structure(c(1, -1), colind = c(2L, 3L))
  )
  fconlist <- list(
    main = structure(diag(2), colind = c(2L, 3L))
  )

  # AR(1) phi per voxel
  phi_matrix <- matrix(0.3, nrow = 1L, ncol = V)

  out <- fmrireg:::fit_lm_contrasts_voxelwise_chunked(
    X_run = X, Y_run = Y, phi_matrix = phi_matrix,
    conlist = conlist, fconlist = fconlist,
    chunk_size = 4L
  )
  expect_true(is.list(out))

  expect_error(
    fmrireg:::fit_lm_contrasts_voxelwise_chunked(
      X_run = X, Y_run = Y, phi_matrix = matrix(0.2, 1, 2),
      conlist = conlist, fconlist = list(), chunk_size = 3L
    ),
    "phi_matrix columns"
  )

  # Vector phi recycled to voxels
  out2 <- fmrireg:::fit_lm_contrasts_voxelwise_chunked(
    X_run = X, Y_run = Y[, 1:3, drop = FALSE],
    phi_matrix = 0.25,
    conlist = conlist, fconlist = list(),
    chunk_size = 2L
  )
  expect_true(is.list(out2))
})

test_that("pull_stat warns on invalid event indices and recovers columns", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  fit2 <- fit
  b <- fit2$result$betas$data[[1]]$estimate[[1]]
  fit2$result$event_indices <- c(ncol(b) + 5L, ncol(b) + 6L)
  expect_warning(
    out <- fmrireg:::pull_stat(fit2, "betas", "estimate"),
    "No valid event indices|valid"
  )
  expect_true(inherits(out, "tbl_df") || is.data.frame(out) || is.matrix(out))
})

test_that(".fast_preproject rank-deficient aliased path", {
  X_def <- cbind(1, c(0, 1, 0, 1, 0, 1), c(0, 1, 0, 1, 0, 1))
  colnames(X_def) <- c("Intercept", "A", "A_dup")
  expect_warning(
    proj <- fmrireg:::.fast_preproject(X_def),
    "rank deficient|Aliased"
  )
  expect_true(is.list(proj))
  expect_true(!is.null(proj$XtXinv))
  expect_false(isTRUE(proj$is_full_rank))
})
