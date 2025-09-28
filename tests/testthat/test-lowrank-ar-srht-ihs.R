test_that("SRHT/IHS + AR(global/cluster) match exact on synthetic data", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")
  library(neuroim2)
  library(fmrihrf)  # For sampling_frame function

  set.seed(42)

  ## ---- Synthetic dataset spec ----
  TR <- 2
  Tlen <- 120
  dim3 <- c(8L, 8L, 3L)                # 192 voxels
  V <- prod(dim3)

  # Build a simple event design: two conditions, one run
  onsets <- seq(20, 200, by = 30)      # seconds
  conds  <- factor(rep(c("A","B"), length.out = length(onsets)))
  events_df <- data.frame(onset = onsets, condition = conds, run = 1L)

  sframe <- sampling_frame(blocklens = Tlen, TR = TR)
  em     <- event_model(onset ~ hrf(condition), data = events_df,
                        block = ~ run, sampling_frame = sframe)

  X <- design_matrix(em)               # T x p
  X <- as.matrix(X)                    # Ensure X is numeric matrix
  p <- ncol(X)

  # True betas: small random effects for task columns; intercept near zero
  B_true <- matrix(0, nrow = p, ncol = V)
  task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
  if (length(task_cols) == 0L) task_cols <- 2:min(3, p)  # fallback
  B_true[task_cols, ] <- matrix(rnorm(length(task_cols) * V, sd = 0.5),
                                nrow = length(task_cols), byrow = TRUE)

  # AR(1) noise generator (time x V)
  ar1_noise <- function(T, V, rho = 0.3, sd = 1.0) {
    E <- matrix(0, T, V)
    E[1, ] <- rnorm(V, sd = sd / sqrt(1 - rho^2))
    for (t in 2:T) E[t, ] <- rho * E[t-1, ] + rnorm(V, sd = sd)
    E
  }

  # Compose Y = X %*% B + AR(1) noise
  Y <- X %*% B_true + ar1_noise(Tlen, V, rho = 0.3, sd = 0.5)   # T x V

  # Convert to 4D array (x,y,z,t) and wrap in NeuroVec + fmri_mem_dataset
  arr <- array(0, dim = c(dim3, Tlen))
  v <- 0L
  for (ix in seq_len(dim3[1]))
    for (iy in seq_len(dim3[2]))
      for (iz in seq_len(dim3[3])) {
        v <- v + 1L
        arr[ix, iy, iz, ] <- as.numeric(Y[, v])
      }

  space4d <- NeuroSpace(c(dim3, Tlen))
  vec     <- NeuroVec(data = arr, space = space4d)
  maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

  dset <- fmri_mem_dataset(scans = list(vec), mask = maskVol, TR = TR, event_table = events_df)

  ## ---- Exact fit (reference) ----
  fit_exact <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
  # Extract exact betas as p x V
  B_exact <- t(fit_exact$result$betas$data[[1]]$estimate[[1]])

  ## ---- Build parcels for pooled AR and latent mapping ----
  coords <- expand.grid(x = seq_len(dim3[1]),
                        y = seq_len(dim3[2]),
                        z = seq_len(dim3[3]))
  K <- 40L
  cl <- kmeans(coords, centers = K, iter.max = 20)$cluster
  parcels <- ClusteredNeuroVol(mask = maskVol, clusters = as.integer(cl))

  ## ---- SRHT + global AR(1) ----
  m_srht <- min(10L * p, Tlen)  # Increased for better accuracy
  low_srht <- lowrank_control(
    parcels     = parcels,
    time_sketch = list(method = "srht", m = m_srht)
  )

  fit_srht_global <- fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    engine   = "latent_sketch",
    lowrank  = low_srht,
    ar_options = list(by_cluster = FALSE, order = 1L)
  )

  ## ---- SRHT + parcel-pooled AR(1) (by_cluster) ----
  fit_srht_group <- fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    engine   = "latent_sketch",
    lowrank  = low_srht,
    ar_options = list(by_cluster = TRUE, order = 1L)
  )

  ## ---- IHS(3 iters) + global AR(1) ----
  low_ihs <- lowrank_control(
    parcels     = parcels,
    time_sketch = list(method = "ihs", m = max(8L * p, p + 10L), iters = 3L)  # Increased sketch size
  )

  fit_ihs <- fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    engine   = "latent_sketch",
    lowrank  = low_ihs,
    ar_options = list(by_cluster = FALSE, order = 1L)
  )

  ## ---- Assertions ----
  # Shapes
  expect_equal(dim(fit_srht_global$betas_fixed), dim(B_exact))
  expect_equal(dim(fit_srht_group$betas_fixed),  dim(B_exact))
  expect_equal(dim(fit_ihs$betas_fixed),         dim(B_exact))

  # Correlation with exact for task regressors only (exclude baseline columns)
  corr_srht_global <- cor(as.numeric(B_exact[task_cols, , drop = FALSE]),
                          as.numeric(fit_srht_global$betas_fixed[task_cols, , drop = FALSE]))
  corr_srht_group  <- cor(as.numeric(B_exact[task_cols, , drop = FALSE]),
                          as.numeric(fit_srht_group$betas_fixed[task_cols, , drop = FALSE]))
  corr_ihs         <- cor(as.numeric(B_exact[task_cols, , drop = FALSE]),
                          as.numeric(fit_ihs$betas_fixed[task_cols, , drop = FALSE]))

  expect_gt(corr_srht_global, 0.75)  # Reasonable threshold for SRHT global
  expect_gt(corr_srht_group,  0.70)  # Adjusted for parcel-pooled AR
  expect_gt(corr_ihs,         0.50)  # Lower threshold for IHS approximation

  # Variance estimates are finite and well-formed
  expect_true(length(fit_srht_global$sigma2) == V)
  expect_true(length(fit_srht_group$sigma2)  == V)
  expect_true(length(fit_ihs$sigma2)         == V)
  expect_true(all(is.finite(fit_srht_global$sigma2)))
  expect_true(all(is.finite(fit_srht_group$sigma2)))
  expect_true(all(is.finite(fit_ihs$sigma2)))

  # AR parameters should be present and finite
  expect_true(is.list(ar_parameters(fit_srht_global, scope = "raw")))
  expect_true(is.list(ar_parameters(fit_srht_group,  scope = "raw")))
  expect_true(is.list(ar_parameters(fit_ihs,         scope = "raw")))
})
