context("cfals wrapper")

library(fmrireg)

simulate_cfals_wrapper_data <- function(hrf_basis, noise_sd = 0.05, signal_scale = 1) {
  sf <- sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 15, 30, 45),
    condition = factor(c("A", "A", "B", "B")),
    block = 1
  )
  emod <- event_model(onset ~ hrf(condition), data = events,
                      block = ~ block, sampling_frame = sf)
  reg_lists <- lapply(emod$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sf,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  sample_times <- samples(sf, global = TRUE)
  X_list <- lapply(regs, function(r)
    evaluate(r, sample_times, precision = sf$precision))
  d <- nbasis(hrf_basis)
  k <- length(X_list)
  v <- 2
  h_true <- matrix(rnorm(d * v), d, v) * signal_scale
  beta_true <- matrix(rnorm(k * v), k, v) * signal_scale
  Y <- matrix(0, nrow(sample_times), v)
  for (c in seq_along(X_list)) {
    Y <- Y + (X_list[[c]] %*% h_true) *
      matrix(rep(beta_true[c, ], each = nrow(Y)), nrow(Y), v)
  }
  Y <- Y + matrix(rnorm(length(Y), sd = noise_sd), nrow(Y), v)
  attr(Y, "sampling_frame") <- sf
  list(Y = Y, event_model = emod, X_list = X_list,
       h_true = h_true, beta_true = beta_true, sframe = sf)
}


test_that("fmrireg_hrf_cfals works across HRF bases", {
  bases <- list(HRF_SPMG3, hrfspline_generator(nbasis = 4))
  for (b in bases) {
    dat <- simulate_cfals_wrapper_data(b)
    fit <- fmrireg_hrf_cfals(dat$Y, dat$event_model, b,
                             lam_beta = 0.1, lam_h = 0.1)
    expect_equal(dim(fit$h_coeffs), c(nbasis(b), ncol(dat$Y)))
    expect_equal(dim(fit$beta_amps), c(length(dat$X_list), ncol(dat$Y)))
    recon <- reconstruction_matrix(b, dat$sframe) %*% fit$h_coeffs
    expect_true(all(is.finite(recon)))
  }
})

test_that("fmrireg_cfals wrapper supports multiple methods", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  methods <- c("ls_svd_only", "ls_svd_1als", "cf_als")
  for (m in methods) {
    fit <- fmrireg_cfals(dat$Y, dat$event_model, HRF_SPMG3,
                         method = m, lambda_b = 0.1, lambda_h = 0.1,
                         lambda_init = 0.5, max_alt = 1)
    expect_equal(dim(fit$h_coeffs), c(nbasis(HRF_SPMG3), ncol(dat$Y)))
    expect_equal(dim(fit$beta_amps), c(length(dat$X_list), ncol(dat$Y)))
  }
})

test_that("cfals handles low-signal data", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3, noise_sd = 0.5, signal_scale = 0.01)
  fit <- fmrireg_hrf_cfals(dat$Y, dat$event_model, HRF_SPMG3)
  expect_lt(mean(fit$gof_per_voxel), 0.2)
})

simple_cfals_data_noise <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 3
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v) + matrix(rnorm(n * v, sd = 0.01), n, v)
  list(X_list = X_list, Y = Y, Xbig = Xbig)
}

test_that("cf_als_engine predictions match canonical GLM", {
  dat <- simple_cfals_data_noise()
  res <- cf_als_engine(dat$X_list, dat$Y,
                       lambda_b = 0,
                       lambda_h = 0,
                       fullXtX_flag = FALSE,
                       max_alt = 1)
  n <- nrow(dat$Y)
  v <- ncol(dat$Y)
  pred_cfals <- matrix(0, n, v)
  for (c in seq_along(dat$X_list)) {
    pred_cfals <- pred_cfals + (dat$X_list[[c]] %*% res$h) *
      matrix(rep(res$beta[c, ], each = n), n, v)
  }
  gamma_hat <- chol2inv(chol(crossprod(dat$Xbig))) %*% crossprod(dat$Xbig, dat$Y)
  pred_glm <- dat$Xbig %*% gamma_hat
  expect_equal(pred_cfals, pred_glm, tolerance = 1e-5)
})

test_that("fullXtX argument is forwarded through fmrireg_cfals", {
  dat <- simulate_cfals_wrapper_data(HRF_SPMG3)
  design <- create_fmri_design(dat$event_model, HRF_SPMG3)
  proj <- project_confounds(dat$Y, design$X_list, NULL)
  direct <- ls_svd_1als_engine(proj$X_list, proj$Y,
                               lambda_init = 0,
                               lambda_b = 0.1,
                               lambda_h = 0.1,
                               fullXtX_flag = TRUE,
                               h_ref_shape_norm = design$h_ref_shape_norm)
  wrap <- fmrireg_cfals(dat$Y, dat$event_model, HRF_SPMG3,
                        method = "ls_svd_1als",
                        fullXtX = TRUE,
                        lambda_init = 0,
                        lambda_b = 0.1,
                        lambda_h = 0.1)
  expect_equal(wrap$h_coeffs, direct$h)
  expect_equal(wrap$beta_amps, direct$beta)
})


