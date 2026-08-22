# tests/testthat/test-landmarks-srht.R

test_that("landmark extension preserves a smooth task field with bounded scale error", {
  skip_on_cran()
  library(neuroim2)

  set.seed(123)
  TR <- 2; Tlen <- 100
  dim3 <- c(6L, 6L, 3L)  # 108 voxels
  space4d <- NeuroSpace(c(dim3, Tlen))
  maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

  # Two-condition design
  onsets <- seq(16, 180, by = 24)
  conds  <- factor(rep(c("A","B"), length.out = length(onsets)))
  events_df <- data.frame(onset = onsets, condition = conds, run = 1L)
  sframe <- sampling_frame(blocklens = Tlen, TR = TR)
  em <- event_model(onset ~ hrf(condition), data = events_df, block = ~ run, sampling_frame = sframe)
  X <- design_matrix(em)
  X <- as.matrix(X)  # Ensure X is numeric matrix
  p <- ncol(X)

  # Simulate Y = X B + AR(1) noise
  V <- prod(dim3)
  task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
  B_true <- matrix(0, p, V)
  coords_sim <- expand.grid(
    z = seq_len(dim3[3]), y = seq_len(dim3[2]), x = seq_len(dim3[1])
  )[, c("x", "y", "z")]
  xyz01 <- sweep(as.matrix(coords_sim), 2, 1, "-")
  xyz01 <- sweep(xyz01, 2, pmax(dim3 - 1, 1), "/")
  smooth_maps <- rbind(
    sin(pi * xyz01[, 1]) * cos(pi * xyz01[, 2]) * (0.5 + xyz01[, 3]),
    cos(pi * xyz01[, 1]) * sin(pi * xyz01[, 2]) * (1 - 0.5 * xyz01[, 3])
  )
  for (j in seq_along(task_cols)) {
    B_true[task_cols[j], ] <- smooth_maps[1L + (j - 1L) %% nrow(smooth_maps), ]
  }
  ar1_noise <- function(T, V, rho = 0.4, sd = 0.5) {
    E <- matrix(0, T, V)
    E[1, ] <- rnorm(V, sd = sd/sqrt(1 - rho^2))
    for (t in 2:T) E[t, ] <- rho * E[t-1, ] + rnorm(V, sd = sd)
    E
  }
  Y <- X %*% B_true + ar1_noise(Tlen, V, sd = 0.15)

  arr <- array(0, dim = c(dim3, Tlen))
  v <- 0L
  for (ix in seq_len(dim3[1])) for (iy in seq_len(dim3[2])) for (iz in seq_len(dim3[3])) {
    v <- v + 1L; arr[ix, iy, iz, ] <- as.numeric(Y[, v])
  }
  vec <- NeuroVec(arr, space4d)
  dset <- fmri_mem_dataset(scans = list(vec), mask = maskVol, TR = TR, event_table = events_df)

  # Parcels via an explicitly convergent k-means route.
  coords <- expand.grid(x = seq_len(dim3[1]), y = seq_len(dim3[2]), z = seq_len(dim3[3]))
  parcel_fit <- expect_no_warning(kmeans(
    coords, centers = 30, iter.max = 1000L, nstart = 10L,
    algorithm = "Lloyd"
  ))
  parcels <- ClusteredNeuroVol(maskVol, parcel_fit$cluster)

  # Compare a full-voxel and landmark fit with the same sketch and AR model.
  full_low <- lowrank_control(
    parcels = parcels,
    time_sketch = list(method = "srht", m = min(10L * p, Tlen))
  )
  control <- fmri_lm_control(
    noise = noise_spec("ar1", pooling = "global")
  )
  set.seed(321)
  fit_full <- expect_no_warning(fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    engine = "latent_sketch", lowrank = full_low, control = control
  ))

  # Design-corrected pooled AR preserves more of the voxel-specific temporal
  # structure than the former mean-series shortcut. Use enough spatial support
  # that this remains a test of landmark extension, not an undersampling edge.
  L <- 40L
  low <- lowrank_control(parcels = parcels, landmarks = L, k_neighbors = 8L,
                         time_sketch = list(method = "srht", m = min(10L * p, Tlen)))
  set.seed(321)
  fit_lm <- expect_no_warning(fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    engine = "latent_sketch", lowrank = low, control = control
  ))

  # Checks
  expect_equal(dim(fit_lm$betas_fixed), dim(fit_full$betas_fixed))
  reference <- as.numeric(fit_full$betas_fixed[task_cols, , drop = FALSE])
  landmark <- as.numeric(fit_lm$betas_fixed[task_cols, , drop = FALSE])
  corr <- cor(reference, landmark)
  rmse <- sqrt(mean((reference - landmark)^2))
  nrmse <- rmse / stats::sd(reference)
  expect_gt(corr, 0.85)
  expect_lt(nrmse, 0.60)
  expect_true(length(fit_lm$sigma2) == V)
  expect_true(all(is.finite(fit_lm$sigma2)))
})
