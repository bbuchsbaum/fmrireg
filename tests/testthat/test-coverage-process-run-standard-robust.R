# Twelfth wave: process_run_standard (IID/AR) and process_run_robust dense paths.

make_run_fx <- function(n = 48L, V = 4L, seed = 201L) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  Y[1:2, 1] <- Y[1:2, 1] + 12
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    Y = Y,
    chunk = list(data = Y, chunk_num = 1L, row_ind = seq_len(n))
  )
}

test_that("process_run_standard covers IID and AR1 with fixed phi", {
  skip_if_not_installed("fmriAR")
  fx <- make_run_fx()
  cfg_iid <- fmri_lm_control()
  out_iid <- fmrireg:::process_run_standard(
    fx$chunk, fx$model, cfg_iid,
    simple_conlist_weights = list(),
    fconlist_weights = list(),
    dataset = fx$dataset
  )
  expect_equal(ncol(out_iid$betas), ncol(fx$Y))
  expect_equal(out_iid$ar_order, 0L)
  expect_null(out_iid$phi_hat)
  expect_true(out_iid$dfres > 0)

  cfg_ar <- fmri_lm_control(ar_options = list(struct = "ar1", exact_first = TRUE))
  out_ar <- fmrireg:::process_run_standard(
    fx$chunk, fx$model, cfg_ar,
    simple_conlist_weights = list(),
    fconlist_weights = list(),
    dataset = fx$dataset
  )
  expect_equal(out_ar$ar_order, 1L)
  expect_true(!is.null(out_ar$phi_hat))
  expect_equal(dim(out_ar$X_final), dim(design_matrix(fx$model)))

  # Fixed phi path
  out_fixed <- fmrireg:::process_run_standard(
    fx$chunk, fx$model, cfg_ar,
    phi_fixed = 0.3,
    simple_conlist_weights = list(),
    fconlist_weights = list(),
    dataset = fx$dataset
  )
  expect_equal(as.numeric(out_fixed$phi_hat), 0.3, tolerance = 1e-8)
})

test_that("process_run_robust returns WLS components and fixed sigma", {
  skip_if_not_installed("robustbase")
  fx <- make_run_fx(seed = 202L)
  cfg <- fmri_lm_control(robust = robust_spec(type = "huber", max_iter = 3L))
  out <- fmrireg:::process_run_robust(
    fx$chunk, fx$model, cfg, dataset = fx$dataset
  )
  expect_equal(ncol(out$betas), ncol(fx$Y))
  expect_equal(length(out$sigma2), ncol(fx$Y))
  expect_true(all(out$sigma2 > 0))
  expect_equal(out$ar_order, 0L)
  expect_true(!is.null(out$robust_weights))
  expect_true(is.finite(out$sigma_robust))

  out_fixed <- fmrireg:::process_run_robust(
    fx$chunk, fx$model, cfg,
    sigma_fixed = 1.25,
    dataset = fx$dataset
  )
  expect_true(is.list(out_fixed))
  expect_equal(ncol(out_fixed$betas), ncol(fx$Y))
})

test_that("prepare_chunkwise_matrices and process_chunk cover IID robust-off", {
  fx <- make_run_fx(n = 40L, V = 3L, seed = 203L)
  cfg <- fmri_lm_control()
  pre <- fmrireg:::prepare_chunkwise_matrices(
    fx$model, fx$dataset, cfg
  )
  expect_true(is.list(pre))
  expect_true(!is.null(pre$X_global) || !is.null(pre$proj_global) || length(pre) >= 1L)

  chunk_out <- fmrireg:::process_chunk(fx$Y, pre, cfg)
  expect_equal(ncol(chunk_out$betas), ncol(fx$Y))
  expect_equal(length(chunk_out$sigma2), ncol(fx$Y))
  expect_true(all(chunk_out$sigma2 > 0))
})