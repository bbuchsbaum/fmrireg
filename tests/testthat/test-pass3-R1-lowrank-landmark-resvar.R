# Pass 3, R1-5: the landmark (Nystrom) extension of the latent_sketch engine
# stored sigma2 = W^2 %*% sigma2_landmarks as the voxel resvar/sigma/rss.
# With rows of W summing to one, sum_l W_vl^2 < 1, so this understated the
# voxel residual variance (by ~1/k for k equal neighbours) and inflated every
# t statistic. The voxel data are available on this path, so each voxel's
# residual variance is now computed from its own data -- the quantity a
# dense fit with the same sketch reports.

.p3l_fixture <- function(dims = c(5L, 5L, 2L), n_time = 80L, seed = 51L) {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RANN")
  set.seed(seed)
  V <- prod(dims)
  arr <- array(rnorm(V * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(onset = seq(5, n_time - 10, by = 12),
                     condition = factor(rep(c("A", "B"), length.out = 6)),
                     run = 1L)
  etab <- etab[seq_len(6), ]
  dset <- neurovec_frame(scans = list(scan), mask = mask, TR = 1, event_table = etab)
  emod <- event_model(onset ~ hrf(condition), data = etab, block = ~ run,
                      sampling_frame = frame_sframe(dset))
  bmod <- baseline_model(basis = "poly", degree = 2, sframe = frame_sframe(dset))
  list(model = fmri_model(emod, bmod, dset), dataset = dset, V = V)
}

.p3l_fit <- function(fx, method, landmarks = NULL, cfg = fmri_lm_control(),
                     seed = 1L) {
  low <- list(time_sketch = list(method = method, m = 40L))
  if (!is.null(landmarks)) {
    low <- c(low, list(landmarks = landmarks, k_neighbors = 4L,
                       kmeans_iter_max = 20L, kmeans_nstart = 1L))
  }
  set.seed(seed)   # same sketch operator for the landmark and dense fits
  suppressMessages(fmrireg:::.run_lowrank_engine(fx$model, fx$dataset,
                                                 lowrank = low, cfg = cfg))
}

test_that("landmark fits report each voxel's own residual variance (= dense fit)", {
  fx <- .p3l_fixture()
  # The landmark extension lives on the AR (whitened) path of the engine.
  for (cfg in list(fmri_lm_control(noise = noise_spec("ar1")),
                   fmri_lm_control(noise = noise_spec("ar2", pooling = "global")))) {
    for (method in c("ihs", "gaussian")) {
      den <- .p3l_fit(fx, method, cfg = cfg)
      lmk <- .p3l_fit(fx, method, landmarks = 6L, cfg = cfg)
      # Betas differ (interpolated), the noise scale does not.
      expect_false(isTRUE(all.equal(lmk$betas_fixed, den$betas_fixed)))
      expect_equal(lmk$result$resvar, den$result$resvar, tolerance = 1e-8)
      expect_equal(lmk$result$sigma, den$result$sigma, tolerance = 1e-8)
      expect_equal(lmk$result$rss, den$result$rss, tolerance = 1e-8)
      expect_equal(lmk$sigma2, den$sigma2, tolerance = 1e-8)
      # Standard errors use that scale with the shared coefficient covariance.
      expect_equal(as.matrix(standard_error(lmk)), as.matrix(standard_error(den)),
                   tolerance = 1e-8)
    }
  }
})

test_that("landmark exact (ihs) residual variance equals the exact voxelwise OLS value", {
  fx <- .p3l_fixture(seed = 52L)
  lmk <- .p3l_fit(fx, "ihs", landmarks = 5L,
                  cfg = fmri_lm_control(noise = noise_spec("ar1")))
  X <- as.matrix(design_matrix(fx$model))
  Y <- as.matrix(fmrireg:::.dset_data_matrix(fx$dataset))
  wh <- fmrireg:::ar_whiten_transform(X, Y, unlist(lmk$ar_coef), FALSE)
  X <- wh$X
  rss <- colSums(qr.resid(qr(X), wh$Y)^2)
  expect_equal(lmk$result$resvar, rss / (nrow(X) - qr(X)$rank), tolerance = 1e-8)
  # Unit-variance noise: the voxel residual variance averages ~1 (SE ~0.03);
  # the propagated W^2 sigma^2 averaged well below that.
  expect_equal(mean(lmk$result$resvar), 1, tolerance = 0.1)
})
