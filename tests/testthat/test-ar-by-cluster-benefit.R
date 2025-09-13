test_that("Parcel AR (by_cluster) reduces residual variance vs global AR", {
  skip_on_cran()
  library(neuroim2)

  set.seed(22)
  TR <- 2; Tlen <- 120
  dim3 <- c(6L, 6L, 2L)  # 72 voxels
  space4d <- NeuroSpace(c(dim3, Tlen))
  maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

  # Simple design
  onsets <- seq(20, 200, by = 30)
  conds  <- factor(rep(c("A","B"), length.out = length(onsets)))
  events_df <- data.frame(onset = onsets, condition = conds, run = 1L)
  sframe <- sampling_frame(blocklens = Tlen, TR = TR)
  em <- event_model(onset ~ hrf(condition), data = events_df, block = ~run, sampling_frame = sframe)
  X <- design_matrix(em); X <- as.matrix(X); p <- ncol(X)

  # Two parcels with different AR(1)
  coords <- expand.grid(x = seq_len(dim3[1]), y = seq_len(dim3[2]), z = seq_len(dim3[3]))
  V <- nrow(coords)
  g <- as.integer(coords$x <= median(coords$x)) + 1L
  parcels <- ClusteredNeuroVol(maskVol, g)

  # Simulate Y = X B + AR noise differing by parcel
  task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
  B_true <- matrix(0, p, V); B_true[task_cols, ] <- matrix(rnorm(length(task_cols) * V, sd = 0.5), length(task_cols), byrow = TRUE)
  ar1_noise_g <- function(T, V, grp, rho1, rho2, sd = 0.5) {
    E <- matrix(0, T, V)
    idx1 <- which(grp == 1L); idx2 <- which(grp == 2L)
    # init
    E[1, idx1] <- rnorm(length(idx1), sd = sd / sqrt(1 - rho1^2))
    E[1, idx2] <- rnorm(length(idx2), sd = sd / sqrt(1 - rho2^2))
    for (t in 2:T) {
      E[t, idx1] <- rho1 * E[t-1, idx1] + rnorm(length(idx1), sd = sd)
      E[t, idx2] <- rho2 * E[t-1, idx2] + rnorm(length(idx2), sd = sd)
    }
    E
  }
  Y <- X %*% B_true + ar1_noise_g(Tlen, V, g, rho1 = 0.7, rho2 = 0.2)

  arr <- array(0, dim = c(dim3, Tlen))
  v <- 0L
  for (ix in seq_len(dim3[1])) for (iy in seq_len(dim3[2])) for (iz in seq_len(dim3[3])) {
    v <- v + 1L; arr[ix, iy, iz, ] <- as.numeric(Y[, v])
  }
  vec <- NeuroVec(arr, space4d)
  dset <- fmri_mem_dataset(scans = list(vec), mask = maskVol, TR = TR, event_table = events_df)

  # Parcels + SRHT sketch setup
  low <- lowrank_control(parcels = parcels, time_sketch = list(method = "srht", m = min(8L * p, Tlen)))

  fit_global <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset, engine = "latent_sketch",
                        lowrank = low, ar_options = list(by_cluster = FALSE, order = 1L))
  fit_group  <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset, engine = "latent_sketch",
                        lowrank = low, ar_options = list(by_cluster = TRUE, order = 1L, shrink_c0 = 100L))

  expect_equal(length(fit_global$sigma2), V)
  expect_equal(length(fit_group$sigma2), V)
  # By-cluster should reduce residual variance for this synthetic
  expect_lt(mean(fit_group$sigma2), mean(fit_global$sigma2))
})
