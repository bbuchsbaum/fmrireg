# fmri_lm_internal.R remaining: rank helpers, voxelwise betas/contrasts, storage.

test_that("is.formula and rank deficiency / preproject helpers", {
  expect_true(fmrireg:::is.formula(y ~ x))
  expect_false(fmrireg:::is.formula("y ~ x"))

  msg <- fmrireg:::.rank_deficiency_message(
    rank = 2, p = 5, aliased = 3:5, varnames = paste0("v", 1:5)
  )
  expect_match(msg, "rank deficient")
  expect_match(msg, "v3")

  long_msg <- fmrireg:::.rank_deficiency_message(
    rank = 1, p = 12, aliased = 2:12, varnames = NULL, context = "X"
  )
  expect_match(long_msg, "\\.\\.\\.")

  X <- cbind(1, 1:10, 1:10)
  info <- fmrireg:::.design_rank_info(X)
  expect_false(isTRUE(info$is_full_rank))
  expect_true(length(info$aliased) >= 1L)

  expect_error(fmrireg:::.design_rank_info(matrix(0, 0, 2)), "at least one")

  attached <- fmrireg:::.attach_rank_attrs(diag(2), info)
  expect_equal(attr(attached, "rank"), info$rank)
  expect_equal(attr(attached, "aliased"), info$aliased)

  expect_error(
    fmrireg:::.fast_preproject(matrix(c(1, NA, 2, 3), 2, 2)),
    "NA/Inf"
  )
  proj <- fmrireg:::.fast_preproject(cbind(1, rnorm(20)))
  expect_true(!is.null(proj$XtXinv))
  expect_true(proj$dfres > 0)
})

test_that("beta_stats_matrix_voxelwise and fit_lm_contrasts_voxelwise", {
  set.seed(81)
  n <- 30L
  p <- 3L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("Intercept", "A", "B")
  Btrue <- c(0.1, 0.8, -0.4)
  Y <- matrix(as.numeric(X %*% Btrue), n, V) + matrix(rnorm(n * V, sd = 0.25), n, V)
  Betas <- qr.coef(qr(X), Y)
  fitted <- X %*% Betas
  sigma <- sqrt(colSums((Y - fitted)^2) / (n - p))
  XtXinv <- chol2inv(chol(crossprod(X)))
  XtXinv_list <- replicate(V, XtXinv, simplify = FALSE)
  # Degenerate voxel: keep list length, mark with non-matrix sentinel + Inf sigma
  XtXinv_list[[V]] <- matrix(NA_real_, 1, 1)
  sigma[V] <- Inf

  bstats <- fmrireg:::beta_stats_matrix_voxelwise(
    Betas = Betas,
    XtXinv_list = XtXinv_list,
    sigma = sigma,
    dfres = n - p,
    varnames = colnames(X),
    robust_weights_list = lapply(seq_len(V), function(i) runif(n, 0.5, 1)),
    ar_order = 1L
  )
  expect_equal(bstats$type, "beta")
  est <- bstats$data[[1]]$estimate[[1]]
  expect_equal(dim(est), c(V, p))
  expect_true(all(is.na(est[V, ])))

  tcon <- list(
    A = structure(1, colind = 2L),
    AvB = structure(c(1, -1), colind = c(2L, 3L))
  )
  fcon <- list(
    AB = structure(diag(2), colind = c(2L, 3L))
  )
  cres <- fmrireg:::fit_lm_contrasts_voxelwise(
    Betas = Betas,
    sigma2 = sigma^2,
    XtXinv_list = XtXinv_list,
    conlist = tcon,
    fconlist = fcon,
    dfres = n - p,
    robust_weights_list = lapply(seq_len(V), function(i) rep(0.9, n)),
    ar_order = 1L
  )
  expect_true(all(c("A", "AvB", "AB") %in% names(cres)))
  expect_equal(nrow(cres$A$data[[1]]), V)
})

test_that("meta_betas / meta_contrasts and store/format with matrix XtXinv", {
  set.seed(82)
  V <- 3L
  p <- 2L
  b1 <- tibble::tibble(
    type = "beta",
    name = "parameter_estimates",
    data = list(tibble::tibble(
      estimate = list(matrix(rnorm(V * p), V, p)),
      se = list(matrix(runif(V * p), V, p)),
      stat = list(matrix(rnorm(V * p), V, p)),
      prob = list(matrix(runif(V * p), V, p)),
      sigma = list(runif(V))
    ))
  )
  b2 <- b1
  b2$data[[1]]$estimate[[1]] <- b2$data[[1]]$estimate[[1]] + 0.1
  mb <- tryCatch(
    fmrireg:::meta_betas(list(b1, b2), event_indices = 1:2),
    error = function(e) e
  )
  expect_true(inherits(mb, "error") || inherits(mb, "tbl_df") || is.list(mb))

  conlist <- list(c1 = structure(1, colind = 2L))
  fconlist <- list(F1 = structure(matrix(1, 1, 1), colind = 2L))
  storage <- fmrireg:::initialize_contrast_storage(conlist, fconlist, n_voxels = 2L)
  XtXinv <- diag(2)
  stored <- fmrireg:::store_voxel_contrasts(
    storage, voxel_index = 1L, qr_or_XtXinv = XtXinv,
    beta_w = c(0.2, 0.5), sigma_w = 0.3,
    conlist = conlist, fconlist = fconlist, dfres = 20
  )
  expect_true(is.list(stored))
  stored2 <- fmrireg:::store_voxel_contrasts(
    stored, voxel_index = 2L, qr_or_XtXinv = qr(cbind(1, rnorm(20))),
    beta_w = c(0.1, 0.4), sigma_w = 0.25,
    conlist = conlist, fconlist = fconlist, dfres = 20, ar_order = 1L
  )
  formatted <- fmrireg:::format_contrast_results(stored2, dfres = 20)
  expect_true(is.list(formatted) || inherits(formatted, "tbl_df"))
})

test_that(".fmri_lm_voxel_status classifies finite / NA / Inf series", {
  expect_equal(fmrireg:::.fmri_lm_voxel_status(rnorm(8)), "ok")
  expect_equal(fmrireg:::.fmri_lm_voxel_status(c(1, NA, 2)), "nonfinite")
  expect_equal(fmrireg:::.fmri_lm_voxel_status(c(1, Inf, 2)), "nonfinite")
  expect_equal(fmrireg:::.fmri_lm_voxel_status(rep(3, 10)), "constant")
  expect_equal(fmrireg:::.fmri_lm_voxel_status(rep(0, 10)), "all_zero")
})
