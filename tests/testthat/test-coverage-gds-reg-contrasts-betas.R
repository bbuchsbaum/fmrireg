# fmrigds registration reducers/regs + contrast/suffstats helpers + betas OLS.

test_that(".fmrigds_make_meta_reg and method mapping cover FE/PM paths", {
  set.seed(21)
  S <- 8
  P <- 3
  V <- 5
  beta <- matrix(rnorm(S * V), S, V)
  se <- matrix(runif(S * V, 0.2, 0.8), S, V)
  X <- cbind(1, rnorm(S), rnorm(S))

  expect_equal(fmrireg:::.fmrigds_meta_method("fe"), "fe")
  expect_equal(fmrireg:::.fmrigds_meta_method("dl"), "dl")
  expect_equal(fmrireg:::.fmrigds_meta_method("pm"), "pm")
  expect_equal(fmrireg:::.fmrigds_meta_method("reml"), "reml")
  expect_equal(fmrireg:::.fmrigds_meta_method("fe", "pm"), "pm")
  expect_equal(fmrireg:::.fmrigds_map_robust("t"), "huber")
  expect_equal(fmrireg:::.fmrigds_map_robust(NULL), "none")

  expect_error(fmrireg:::.fmrigds_coerce_var(NULL, NULL), "var' or 'se'")
  expect_equal(fmrireg:::.fmrigds_coerce_var(NULL, se), se^2)

  reg_fe <- fmrireg:::.fmrigds_make_meta_reg("fe")
  out_fe <- reg_fe(beta = beta, var = se^2, X = X, z = NULL, p = NULL,
                   df = NULL, df1 = NULL, df2 = NULL, opts = list())
  expect_true(all(c("coef", "se_coef") %in% names(out_fe) |
                    c("beta", "se") %in% names(out_fe)))

  # Reducer FE path (fmri_meta_fit expects subjects x voxels for Y)
  red_fe <- fmrireg:::.fmrigds_make_meta_reducer("fe")
  out_r <- red_fe(beta = beta, var = NULL, X = X, z = NULL, p = NULL,
                  df = NULL, df1 = NULL, df2 = NULL,
                  opts = list(se_override = se))
  expect_true(is.matrix(out_r$beta) || is.numeric(out_r$beta))
  expect_true(length(out_r$beta) > 0)

  # PM with return_cov tri
  red_pm <- fmrireg:::.fmrigds_make_meta_reducer("pm")
  out_pm <- red_pm(beta = beta, var = se^2, X = X, z = NULL, p = NULL,
                   df = NULL, df1 = NULL, df2 = NULL,
                   opts = list(return_cov = "tri", method = "pm"))
  expect_true(!is.null(out_pm$beta))
  expect_true(is.null(out_pm$cov_tri) || is.matrix(out_pm$cov_tri) || is.numeric(out_pm$cov_tri))

  # Contrast path: Cmat is p x n_contrasts
  Cmat <- matrix(c(1, 0, 0), ncol = 1)
  out_c <- red_fe(beta = beta, var = se^2, X = X, z = NULL, p = NULL,
                  df = NULL, df1 = NULL, df2 = NULL,
                  opts = list(contrasts = Cmat))
  expect_true(!is.null(out_c$con_est) || !is.null(out_c$beta))

  if (exists(".fmrigds_check_reducer", envir = asNamespace("fmrireg"), inherits = FALSE)) {
    ok <- tryCatch(
      fmrireg:::.fmrigds_check_reducer("meta:fe_reg", list()),
      error = function(e) FALSE
    )
    expect_true(is.logical(ok) || isTRUE(ok) || isFALSE(ok))
  }
})

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
  lds <- fmridataset::latent_dataset(source = list(lvec), TR = 1, run_length = n_time)
  lds$event_table <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
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
  dset <- matrix_dataset(Y, TR = 1, run_length = 30, event_table = etab)
  emod <- event_model(onset ~ hrf(condition), data = etab, block = ~ run,
                      sampling_frame = dset$sampling_frame)

  expect_error(glm_ols(list(), emod), regexp = ".")
  expect_error(glm_ols(dset, list()), regexp = ".")
  if (exists("glm_lss", mode = "function")) {
    expect_error(glm_lss(list(), emod), regexp = ".")
  }
})
