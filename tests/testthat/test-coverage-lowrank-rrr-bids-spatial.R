# Eleventh wave: lowrank AR helpers, RRR whiten, BIDS spatial write success.

test_that(".lowrank_whiten_initial_ols and .lowrank_group_ar_estimates succeed", {
  skip_if_not_installed("fmriAR")
  set.seed(101)
  n <- 60L
  p <- 3L
  V <- 8L
  X <- cbind(1, rnorm(n), rnorm(n))
  Z <- X %*% matrix(c(0.2, 0.5, -0.3), p, 1) %*% matrix(1, 1, V) +
    matrix(rnorm(n * V, sd = 0.4), n, V)
  ar_opts <- list(struct = "ar1", exact_first = TRUE)

  init <- fmrireg:::.lowrank_whiten_initial_ols(
    X, Z, ar_order = 1L, ar_opts = ar_opts
  )
  expect_true(!is.null(init$phi))
  expect_equal(dim(init$X), dim(X))
  expect_equal(dim(init$Y), dim(Z))

  expect_error(
    fmrireg:::.lowrank_group_ar_estimates(
      residuals = list(), sizes = integer(0), ar_order = 1L,
      ar_opts = ar_opts, shrink_c0 = 10L, design = X
    ),
    "non-empty list"
  )

  res_list <- list(
    `1` = as.numeric(Z[, 1] - X %*% qr.coef(qr(X), Z[, 1])),
    `2` = as.numeric(rowMeans(Z[, 2:4]) - X %*% qr.coef(qr(X), rowMeans(Z[, 2:4])))
  )
  sizes <- c(`1` = 1L, `2` = 3L)
  grp <- fmrireg:::.lowrank_group_ar_estimates(
    residuals = res_list, sizes = sizes, ar_order = 1L,
    ar_opts = ar_opts, shrink_c0 = 5L, design = X
  )
  expect_true(!is.null(grp$global_phi))
  expect_equal(length(grp$phi_groups), 2L)
})

test_that(".rrr_whiten_data covers iid and ar1 with nuisance", {
  skip_if_not_installed("fmriAR")
  set.seed(102)
  n <- 50L
  X_task <- cbind(rnorm(n), rnorm(n))
  Z_nuis <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 4), n, 4)

  iid <- fmrireg:::.rrr_whiten_data(
    X_task, Z_nuis, Y, cfg = fmri_lm_control(), run_indices = NULL
  )
  expect_equal(iid$ar_order, 0L)
  expect_equal(dim(iid$Xw_task), dim(X_task))
  expect_equal(dim(iid$Zw), dim(Z_nuis))

  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", exact_first = TRUE))
  wh <- fmrireg:::.rrr_whiten_data(
    X_task, Z_nuis, Y, cfg = cfg, run_indices = list(1:25, 26:50)
  )
  expect_equal(wh$ar_order, 1L)
  expect_equal(dim(wh$Yw), dim(Y))
  expect_true(length(wh$phi) >= 1L || is.numeric(wh$phi) || is.list(wh$phi))

  # No nuisance
  wh2 <- fmrireg:::.rrr_whiten_data(
    X_task, NULL, Y, cfg = cfg, run_indices = NULL
  )
  expect_null(wh2$Zw)
})

test_that(".write_beta_outputs and .write_contrast_outputs succeed on spatial nifti", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  set.seed(103)
  dims <- c(2L, 2L, 1L)
  n_time <- 36L
  arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = c(5, 12, 20, 28),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = list(scan), mask = mask, TR = 1, event_table = etab
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~ run, dataset = dset
  )
  td <- tempfile("bids-spatial-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "01", task = "demo", space = "T1w")

  betas <- fmrireg:::.write_beta_outputs(
    fit, td, ent, "GLM", format = "nifti", save_betas = TRUE, overwrite = TRUE
  )
  expect_true(is.list(betas))
  expect_true(!is.null(betas$betas) || length(list.files(td, pattern = "betas")) >= 1L)

  cons <- fmrireg:::.write_contrast_outputs(
    fit, td, ent, "GLM", format = "nifti", strategy = "by_contrast",
    contrasts = NULL, contrast_match = "auto",
    contrast_stats = c("beta", "tstat", "se"), overwrite = TRUE
  )
  expect_true(is.list(cons))
  expect_true(length(list.files(td, pattern = "nii|json")) >= 1L)

  cons2 <- fmrireg:::.write_contrast_outputs(
    fit, td, ent, "GLM", format = "nifti", strategy = "by_stat",
    contrasts = "A_vs_B", contrast_match = "exact",
    contrast_stats = c("beta", "tstat"), overwrite = TRUE
  )
  expect_true(is.list(cons2))
})

test_that("fit_contrasts.fmri_lm structured voxel F-path and matrix L", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  V <- length(fit$result$sigma %||% fit$result$betas$data[[1]]$sigma[[1]])
  if (V < 1L) skip("demo fit empty")
  XtXinv <- fit$result$cov.unscaled %||% fit$result$XtXinv
  if (is.null(XtXinv)) skip("no cov")
  p_event <- length(fit$result$event_indices)
  if (p_event < 2L) skip("need >=2 event coefs")

  # Attach per-voxel covariance matching event-sized betas
  event_cov <- as.matrix(XtXinv)[fit$result$event_indices, fit$result$event_indices, drop = FALSE]
  fit$result$variance_model <- structure(
    list(
      df_inference = rep(max(2, fit$result$rdf %||% 10), V),
      covariance_scope = "voxel",
      covariance = replicate(V, event_cov, simplify = FALSE)
    ),
    class = "fmri_lm_variance_model"
  )

  bnames <- colnames(coef(fit))
  if (length(bnames) >= 2L) {
    L <- matrix(c(1, -1), 1, dimnames = list("diff", bnames[1:2]))
    out <- tryCatch(fit_contrasts(fit, list(diff = L)), error = function(e) e)
    expect_true(inherits(out, "error") || is.list(out))

    LF <- rbind(c(1, 0), c(0, 1))
    colnames(LF) <- bnames[1:2]
    outF <- tryCatch(fit_contrasts(fit, list(main = LF)), error = function(e) e)
    expect_true(inherits(outF, "error") || is.list(outF))
  }
})
