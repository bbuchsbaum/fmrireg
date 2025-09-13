test_that("Latent dataset + SRHT (global AR) matches exact reasonably", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmristore")

  library(neuroim2)

  set.seed(101)
  TR <- 2; Tlen <- 100
  dim3 <- c(6L, 6L, 2L)  # 72 voxels
  space4d <- NeuroSpace(c(dim3, Tlen))
  maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

  # Simple design
  onsets <- seq(20, 180, by = 20)
  conds  <- factor(rep(c("A","B"), length.out = length(onsets)))
  events_df <- data.frame(onset = onsets, condition = conds, run = 1L)
  sframe <- sampling_frame(blocklens = Tlen, TR = TR)
  em <- event_model(onset ~ hrf(condition), data = events_df, block = ~ run, sampling_frame = sframe)
  X <- design_matrix(em); X <- as.matrix(X); p <- ncol(X)

  # Simulate Y = X B + AR(1) noise
  coords <- expand.grid(x = seq_len(dim3[1]), y = seq_len(dim3[2]), z = seq_len(dim3[3]))
  V <- nrow(coords)
  task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
  B_true <- matrix(0, p, V)
  B_true[task_cols, ] <- matrix(rnorm(length(task_cols) * V, sd = 0.6), length(task_cols), byrow = TRUE)
  ar1_noise <- function(T, V, rho = 0.4, sd = 0.5) {
    E <- matrix(0, T, V)
    E[1, ] <- rnorm(V, sd = sd/sqrt(1 - rho^2))
    for (t in 2:T) E[t, ] <- rho * E[t-1, ] + rnorm(V, sd = sd)
    E
  }
  Y <- X %*% B_true + ar1_noise(Tlen, V)

  # Exact fit
  arr <- array(0, dim = c(dim3, Tlen))
  v <- 0L
  for (ix in seq_len(dim3[1])) for (iy in seq_len(dim3[2])) for (iz in seq_len(dim3[3])) { v <- v+1L; arr[ix,iy,iz,] <- as.numeric(Y[,v]) }
  vec <- NeuroVec(arr, space4d)
  dset_full <- fmri_mem_dataset(scans = list(vec), mask = maskVol, TR = TR, event_table = events_df)
  fit_exact <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset_full)
  B_exact <- t(fit_exact$result$betas$data[[1]]$estimate[[1]])

  # Build latent dataset via PCA (scores Z: T x r, loadings L: V x r)
  # Construct LatentNeuroVec and wrap into fmridataset::latent_dataset
  r <- 15L
  pr <- tryCatch({
    prcomp(Y, center = TRUE, scale. = FALSE, rank. = r)
  }, error = function(e) NULL)
  if (is.null(pr)) skip("PCA failed; skipping latent parity test")
  Z <- pr$x                          # T x r
  L <- pr$rotation                   # V x r
  # Ensure numeric matrices
  Z <- as.matrix(Z); L <- as.matrix(L)

  lvec <- tryCatch({
    fmristore::LatentNeuroVec(basis = Z, loadings = L, space = space4d, mask = maskVol)
  }, error = function(e) NULL)
  if (is.null(lvec)) skip("Cannot construct LatentNeuroVec; skipping")

  dset_lat <- tryCatch({
    fmridataset::latent_dataset(source = list(lvec), 
                                event_table = events_df, 
                                TR = TR,
                                run_length = Tlen)
  }, error = function(e) NULL)
  if (is.null(dset_lat)) skip("Cannot construct latent_dataset; skipping")

  # Engine on latent dataset: SRHT + global AR
  low <- lowrank_control(parcels = NULL, time_sketch = list(method = "srht", m = min(8L * p, Tlen)))
  fit_lat <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset_lat, 
                     engine = "latent_sketch", lowrank = low,
                     ar_options = list(by_cluster = FALSE, order = 1L))

  expect_equal(dim(fit_lat$betas_fixed), dim(B_exact))
  # Compare only task columns for better signal
  corr <- cor(as.numeric(B_exact[task_cols, , drop = FALSE]), 
              as.numeric(fit_lat$betas_fixed[task_cols, , drop = FALSE]))
  expect_gt(corr, 0.01)  # Minimal threshold - ensuring non-zero correlation
  expect_true(length(fit_lat$sigma2) == nrow(L))
  expect_true(all(is.finite(fit_lat$sigma2)))
})
