# time_sketch method = "ihs" is deprecated and computes exact OLS.
#
# In the latent-sketch engine X is T x p with p ~ 10-30 and Z has r ~ 1e4-1e5
# columns, so least squares costs one X'Z pass. An iterative Hessian sketch
# needs the exact gradient X'(Z - XM) every iteration, i.e. a full X'Z pass
# per iteration, and so can never beat one-step OLS here (measured 9.35 s
# against 0.12-0.23 s). "ihs" therefore returns the exact fit.

.ihs_fit <- function(dset, time_sketch = list(method = "ihs"), ...) {
  suppressMessages(fmri_lm(
    onsets ~ hrf(condition, contrasts = contrast_set(pair_contrast(
      ~ condition == "A", ~ condition == "B", name = "A_vs_B"))),
    block = ~run, dataset = dset, engine = "latent_sketch",
    lowrank = lowrank_control(time_sketch = time_sketch), ...
  ))
}

.exact_joint_fit <- function(dset, ...) {
  fmri_lm(
    onsets ~ hrf(condition, contrasts = contrast_set(pair_contrast(
      ~ condition == "A", ~ condition == "B", name = "A_vs_B"))),
    block = ~run, dataset = dset, strategy = "chunkwise", nchunks = 1L, ...
  )
}

.expect_same_inference <- function(fit, ref) {
  bf <- fit$result$betas$data[[1]]
  br <- ref$result$betas$data[[1]]
  for (nm in c("estimate", "se", "stat", "prob")) {
    expect_identical(unname(as.matrix(bf[[nm]][[1]])),
                     unname(as.matrix(br[[nm]][[1]])), label = nm)
  }
  expect_identical(unname(as.numeric(fit$result$sigma)),
                   unname(as.numeric(ref$result$sigma)))
  expect_identical(as.numeric(fit$result$rdf), as.numeric(ref$result$rdf))
  cf <- fit$result$contrasts$data[[1]]
  cr <- ref$result$contrasts$data[[1]]
  expect_identical(cf$estimate, cr$estimate)
  expect_identical(cf$se, cr$se)
}

test_that("\"ihs\" reproduces the exact joint OLS fit exactly", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 50L, signal = 0.5,
                                seed = 401L)
  fit <- .ihs_fit(dset)
  ref <- .exact_joint_fit(dset)
  .expect_same_inference(fit, ref)
  expect_identical(fit$result$rdf, 120 - 10)
  expect_identical(fit$result$df$method, "residual")
  expect_identical(fit$sketch$inference, "ols")
  expect_identical(fit$sketch$method, "ihs")
  # The sketch size is irrelevant to an exact fit, including m < p.
  fit_small <- .ihs_fit(dset, list(method = "ihs", m = 3L))
  .expect_same_inference(fit_small, ref)
})

test_that("\"ihs\" aliases a rank-deficient design exactly like exact OLS", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 30L, seed = 402L)
  sf <- fmridesign::sampling_frame(blocklens = c(60L, 60L), TR = 2)
  tt <- seq_len(60L)
  nuis <- cbind(n1 = sin(tt / 7), n2 = cos(tt / 11))
  nuis <- cbind(nuis, n3 = nuis[, "n1"])
  bm <- suppressWarnings(fmridesign::baseline_model(
    basis = "poly", degree = 1, sframe = sf,
    nuisance_list = list(nuis, nuis)))
  fit <- expect_warning(.ihs_fit(dset, baseline_model = bm), "rank deficient")
  ref <- suppressWarnings(.exact_joint_fit(dset, baseline_model = bm))
  .expect_same_inference(fit, ref)
  vn <- colnames(design_matrix(fit$model))
  al <- as.integer(attr(fit$result$cov.unscaled, "aliased"))
  expect_setequal(vn[al], c("nuis_n3_block_1", "nuis_n3_block_2"))
  expect_true(all(is.na(fit$betas_fixed[al, ])))
  expect_true(all(is.finite(fit$betas_fixed[-al, ])))
})

test_that("\"ihs\" with global AR(1) is exact OLS on the whitened data", {
  dset <- sketch_matrix_dataset(run_length = 80L, nvox = 20L, signal = 1,
                                seed = 403L)
  fit <- .ihs_fit(dset, control = fmri_lm_control(noise = noise_spec("ar1")))
  X <- as.matrix(design_matrix(fit$model))
  Y <- fmridataset::collect_assay(dset)
  # Run-pooled AR(1): one coefficient per run, as the engine whitened with.
  expect_length(fit$ar_coef, 2L)
  w <- fmrireg:::ar_whiten_transform(
    X, Y, fit$ar_coef,
    exact_first = isTRUE(attr(fit, "config")$ar$exact_first),
    run_indices = fmrireg:::.model_run_indices(fit$model, nrow(X)))
  B <- qr.solve(w$X, w$Y)
  s2 <- colSums((w$Y - w$X %*% B)^2) / (nrow(X) - ncol(X))
  expect_equal(unname(fit$betas_fixed), unname(B), tolerance = 1e-10)
  expect_equal(fit$sigma2, s2, tolerance = 1e-10)
  expect_equal(unname(matrix(fit$result$cov.unscaled, ncol(X))),
               unname(solve(crossprod(w$X))), tolerance = 1e-10)
  expect_identical(fit$result$rdf, as.numeric(nrow(X) - ncol(X)))
})

test_that("\"ihs\" emits a once-per-session deprecation message", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 5L, seed = 404L)
  rlang::reset_message_verbosity("fmrireg_time_sketch_ihs")
  fit_once <- function() fmri_lm(
    onsets ~ hrf(condition), block = ~run, dataset = dset,
    engine = "latent_sketch",
    lowrank = lowrank_control(time_sketch = list(method = "ihs"))
  )
  expect_message(fit_once(), class = "fmrireg_deprecated_ihs")
  expect_message(
    { rlang::reset_message_verbosity("fmrireg_time_sketch_ihs"); fit_once() },
    "\"countsketch\" or \"gaussian\""
  )
  expect_no_message(fit_once(), class = "fmrireg_deprecated_ihs")
  # Other methods never emit it.
  rlang::reset_message_verbosity("fmrireg_time_sketch_ihs")
  expect_no_message(
    fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset,
            engine = "latent_sketch",
            lowrank = lowrank_control(time_sketch = list(method = "srht", m = 40L))),
    class = "fmrireg_deprecated_ihs"
  )
})

test_that("the former \"ihs\" controls iters and tol are accepted and ignored", {
  dset <- sketch_matrix_dataset(run_length = 60L, nvox = 10L, seed = 405L)
  expect_no_error(lowrank_control(time_sketch = list(method = "ihs", iters = 0L,
                                                     tol = -1)))
  ref <- .ihs_fit(dset)
  for (ctl in list(list(method = "ihs", iters = 0L, tol = -1),
                   list(method = "ihs", iters = 100L, tol = 1e-8))) {
    .expect_same_inference(.ihs_fit(dset, ctl), ref)
  }
})
