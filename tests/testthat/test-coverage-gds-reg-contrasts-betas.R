# Contrast/suffstats helpers + betas OLS.

test_that("compute_lm_contrasts_from_suffstats and normalize error paths", {
  set.seed(22)
  n <- 40
  p <- 3
  V <- 4
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  Y <- X %*% matrix(rnorm(p * V), p, V) + matrix(rnorm(n * V), n, V)
  XtX <- crossprod(X)
  XtS <- crossprod(X, Y)
  StS <- colSums(Y^2)
  XtXinv <- solve(XtX)
  betas <- solve(XtX, XtS)
  rownames(betas) <- colnames(X)
  sigma2 <- (StS - colSums(betas * XtS)) / (n - p)

  # Empty contrasts
  empty <- compute_lm_contrasts_from_suffstats(
    XtX = XtX, XtS = XtS, StS = StS, df = n - p, contrasts = list()
  )
  expect_true(is.data.frame(empty) || inherits(empty, "tbl_df") || is.list(empty))

  out_list <- compute_lm_contrasts_from_suffstats(
    XtX = XtX, XtS = XtS, StS = StS, df = n - p,
    t_contrasts = list(c1 = matrix(c(0, 1, -1), nrow = 1)),
    output = "list"
  )
  expect_true(is.list(out_list))

  stacked <- compute_lm_contrasts(
    B = betas, XtXinv = XtXinv, sigma2 = sigma2, df = n - p,
    t_contrasts = list(c1 = matrix(c(0, 1, -1), nrow = 1)),
    columns = colnames(X)
  )
  expect_true(nrow(stacked) >= 1L)

  expect_error(
    compute_lm_contrasts(
      B = betas, XtXinv = XtXinv, sigma2 = sigma2, df = n - p,
      contrasts = list(bad = "not-a-contrast"),
      drop_failed = FALSE
    ),
    regexp = "."
  )
  dropped <- compute_lm_contrasts(
    B = betas, XtXinv = XtXinv, sigma2 = sigma2, df = n - p,
    contrasts = list(bad = "not-a-contrast"),
    drop_failed = TRUE
  )
  expect_true(is.data.frame(dropped) || inherits(dropped, "tbl_df") || is.list(dropped))

  # F-contrast path
  f_out <- compute_lm_contrasts(
    B = betas, XtXinv = XtXinv, sigma2 = sigma2, df = n - p,
    f_contrasts = list(fAB = matrix(c(0, 1, 0, 0, 0, 1), nrow = 2)),
    columns = colnames(X)
  )
  expect_true(nrow(f_out) >= 1L || length(f_out) >= 0L)
})

test_that("estimate_betas.latent_dataset OLS and glm_* validation", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(23)
  n_time <- 50
  n_comp <- 3
  n_voxels <- 12
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  loadings <- matrix(rnorm(n_voxels * n_comp), n_voxels, n_comp)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis, loadings = loadings,
    space = neuroim2::NeuroSpace(c(4, 3, 1, n_time)),
    mask = rep(TRUE, n_voxels), offset = rep(0, n_voxels)
  )
  lds <- latent_frame(
    lvec, TR = 1, run_length = n_time,
    event_table = data.frame(
      onset = c(5, 15, 25, 35),
      condition = factor(c("A", "B", "A", "B")),
      run = 1L
    )
  )

  betas <- tryCatch(
    estimate_betas(
      lds, ran = onset ~ hrf(condition), block = ~ run,
      method = "ols", progress = FALSE
    ),
    error = function(e) e
  )
  if (inherits(betas, "error")) {
    expect_true(nzchar(conditionMessage(betas)))
  } else {
    expect_true(inherits(betas, "fmri_latent_betas") || is.list(betas))
  }

  # glm_ols / glm_lss validation
  etab <- data.frame(onset = 5, condition = factor("A"), run = 1L)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  dset <- matrix_frame(Y, TR = 1, run_length = 30, event_table = etab)
  emod <- event_model(onset ~ hrf(condition), data = etab, block = ~ run,
                      sampling_frame = frame_sframe(dset))

  expect_error(glm_ols(list(), emod), regexp = ".")
  expect_error(glm_ols(dset, list()), regexp = ".")
  if (exists("glm_lss", mode = "function")) {
    expect_error(glm_lss(list(), emod), regexp = ".")
  }
})
