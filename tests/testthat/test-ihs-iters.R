test_that("IHS with more iterations improves or matches correlation vs exact", {
  skip_on_cran()
  library(neuroim2)

  set.seed(202)
  TR <- 2; Tlen <- 120
  dim3 <- c(6L, 6L, 3L)  # 108 voxels
  space4d <- NeuroSpace(c(dim3, Tlen))
  maskVol <- LogicalNeuroVol(array(TRUE, dim3), NeuroSpace(dim3))

  # Two-condition design
  onsets <- seq(18, 200, by = 24)
  conds  <- factor(rep(c("A","B"), length.out = length(onsets)))
  events_df <- data.frame(onset = onsets, condition = conds, run = 1L)
  sframe <- sampling_frame(blocklens = Tlen, TR = TR)
  em <- event_model(onset ~ hrf(condition), data = events_df, block = ~ run, sampling_frame = sframe)
  X <- design_matrix(em); X <- as.matrix(X); p <- ncol(X)

  # Simulate Y = X B + AR(1) noise
  V <- prod(dim3)
  task_cols <- which(grepl("condition|hrf", colnames(X), ignore.case = TRUE))
  B_true <- matrix(0, p, V)
  B_true[task_cols, ] <- matrix(rnorm(length(task_cols) * V, sd = 0.5), length(task_cols), byrow = TRUE)
  # Express the effect by its peak height relative to noise, not as a raw
  # coefficient: the HRF's raw scale belongs to fmrihrf (it changed ~10x when
  # the SPMG undershoot was corrected). A peak of 1.75 per unit beta is the
  # scale this fixture's SNR was calibrated against.
  B_true[task_cols, ] <- sweep(B_true[task_cols, , drop = FALSE], 1,
                               1.75 / apply(abs(X[, task_cols, drop = FALSE]), 2, max), `*`)
  ar1_noise <- function(T, V, rho = 0.35, sd = 0.5) {
    E <- matrix(0, T, V)
    E[1, ] <- rnorm(V, sd = sd/sqrt(1 - rho^2))
    for (t in 2:T) E[t, ] <- rho * E[t-1, ] + rnorm(V, sd = sd)
    E
  }
  Y <- X %*% B_true + ar1_noise(Tlen, V)

  arr <- array(0, dim = c(dim3, Tlen))
  v <- 0L
  for (ix in seq_len(dim3[1])) for (iy in seq_len(dim3[2])) for (iz in seq_len(dim3[3])) { v <- v+1L; arr[ix,iy,iz,] <- as.numeric(Y[,v]) }
  vec <- NeuroVec(arr, space4d)
  dset <- neurovec_frame(scans = list(vec), mask = maskVol, TR = TR, event_table = events_df)

  # Exact
  fit_exact <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)
  B_exact <- t(fit_exact$result$betas$data[[1]]$estimate[[1]])

  # IHS iters 1 vs 3
  low1 <- lowrank_control(parcels = NULL, time_sketch = list(method = "ihs", m = max(8L * p, p + 10L), iters = 1L))
  fit1 <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset, engine = "latent_sketch", lowrank = low1,
                  ar_options = list(by_cluster = FALSE, order = 1L))

  low3 <- lowrank_control(parcels = NULL, time_sketch = list(method = "ihs", m = max(8L * p, p + 10L), iters = 3L))
  fit3 <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset, engine = "latent_sketch", lowrank = low3,
                  ar_options = list(by_cluster = FALSE, order = 1L))

  corr1 <- cor(as.numeric(B_exact), as.numeric(fit1$betas_fixed))
  corr3 <- cor(as.numeric(B_exact), as.numeric(fit3$betas_fixed))
  expect_true(corr1 > 0.40)  # Realistic threshold for single IHS iteration
  expect_gte(corr3, corr1)  # no worse, typically better
})

test_that("IHS kernel converges to the exact least-squares solution", {
  # Iterative Hessian sketching sketches only the Hessian; the gradient uses
  # the full data, so iterates converge to the exact OLS solution rather than
  # to (or around) the solution of one sketched problem.
  set.seed(77)
  Tlen <- 120L
  X <- cbind(1, rnorm(Tlen), sin(seq_len(Tlen) / 5))
  Z <- X %*% matrix(rnorm(3 * 5), 3) + matrix(rnorm(Tlen * 5), Tlen)
  B_ols <- qr.solve(X, Z)
  rss <- function(M) sum((Z - X %*% M)^2)

  errs <- vapply(c(1L, 3L, 25L), function(it) {
    set.seed(78)
    max(abs(ihs_latent_solve(X, Z, m = 24L, iters = it)$M - B_ols))
  }, numeric(1))
  expect_lt(errs[3], 1e-6)
  expect_true(all(diff(errs) < 0))

  rss_path <- vapply(1:6, function(it) {
    set.seed(78)
    rss(ihs_latent_solve(X, Z, m = 24L, iters = it)$M)
  }, numeric(1))
  expect_true(all(diff(rss_path) <= 1e-10 * rss_path[-1]))
  expect_lt(rss_path[6], rss_path[1])
})
