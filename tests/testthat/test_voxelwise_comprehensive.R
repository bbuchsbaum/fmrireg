context("Comprehensive voxelwise AR contrasts")

set.seed(123)

# Helper data setup
etab <- data.frame(onset = c(1, 10), cond = factor(c("A", "B")), run = c(1, 1))
con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"))


test_that("single voxel voxelwise matches global AR", {
  y <- as.numeric(arima.sim(model = list(ar = 0.4), n = 20))
  dset <- matrix_dataset(matrix(y, ncol = 1), TR = 1, run_length = 20,
                         event_table = etab)

  mod_global <- fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~ run,
                        dataset = dset, use_fast_path = FALSE,
                        cor_struct = "ar1", ar_voxelwise = FALSE)
  mod_vox <- fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~ run,
                     dataset = dset, use_fast_path = FALSE,
                     cor_struct = "ar1", ar_voxelwise = TRUE)
  expect_equal(coef(mod_vox), coef(mod_global), tolerance = 1e-1)
  expect_equal(stats(mod_vox, "contrasts"), stats(mod_global, "contrasts"),
               tolerance = 5e-1)
})

test_that("chunked voxelwise engine matches direct computation", {
  n <- 30
  p <- 3
  V <- 6
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  phi <- matrix(runif(1 * V, -0.3, 0.3), 1, V)
  Y <- matrix(0, n, V)
  Betas <- matrix(0, p, V)
  sigma <- numeric(V)
  XtXinv_list <- vector("list", V)

  for (v in seq_len(V)) {
    noise <- as.numeric(arima.sim(n = n, model = list(ar = phi[, v]), sd = 0.1))
    Y[, v] <- X %*% rnorm(p) + noise
    tmp <- ar_whiten_transform(X, Y[, v, drop = FALSE], phi[, v], exact_first = TRUE)
    qr_w <- qr(tmp$X)
    Betas[, v] <- qr.coef(qr_w, tmp$Y)
    sigma[v] <- sqrt(sum(qr.resid(qr_w, tmp$Y)^2) / (nrow(tmp$X) - qr_w$rank))
    XtXinv_list[[v]] <- chol2inv(qr.R(qr_w))
  }

  # Use consistent contrast setup: scalar contrast with scalar colind
  clist <- list(A = 1)  # Apply weight 1 to column 2
  attr(clist$A, "colind") <- 2L
  fclist <- list()
  dfres <- nrow(X) - qr(X)$rank

  res_direct <- fit_lm_contrasts_voxelwise(Betas, sigma^2, XtXinv_list,
                                           clist, fclist, dfres)
  res_chunk <- fit_lm_contrasts_voxelwise_chunked(X, Y, phi, clist, fclist,
                                                  chunk_size = 2)
  expect_equal(res_chunk$A$data[[1]], res_direct$A$data[[1]], tolerance = 1e-6)
})

test_that("constant voxel handled without error", {
  n <- 20
  X <- cbind(1, rnorm(n))
  Y <- matrix(1, n, 1)
  phi <- 0.2
  # Use consistent contrast setup: scalar contrast with scalar colind
  clist <- list(A = 1)  # Apply weight 1 to column 2
  attr(clist$A, "colind") <- 2L
  res <- fit_lm_contrasts_voxelwise_chunked(X, Y, phi, clist, list(),
                                            chunk_size = 1)
  stat <- res$A$data[[1]]$stat
  expect_true(is.finite(stat) || is.nan(stat))
})

test_that("effective df adjusts p-values for AR order", {
  n <- 25
  X <- cbind(1, rnorm(n))
  phi <- 0.3
  y <- as.numeric(arima.sim(n = n, model = list(ar = phi), sd = 0.1))
  # Use consistent contrast setup: scalar contrast with scalar colind
  clist <- list(A = 1)  # Apply weight 1 to column 2
  attr(clist$A, "colind") <- 2L
  res <- fit_lm_contrasts_voxelwise_chunked(X, matrix(y, ncol = 1), phi,
                                            clist, list(), chunk_size = 1)
  tval <- res$A$data[[1]]$stat
  dfres <- n - qr(X)$rank
  df_eff <- calculate_effective_df(n, ncol(X), ar_order = 1)
  p_expected <- 2 * pt(-abs(tval), df_eff)
  expect_equal(res$A$data[[1]]$prob, p_expected, tolerance = 1e-6)
})
