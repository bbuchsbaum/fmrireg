test_that("process_run_integrated honors sigma_fixed in robust pipeline", {
  set.seed(1001)
  n <- 120
  p <- 4
  v <- 20
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rt(n * v, df = 3), n, v)

  cfg <- fmri_lm_control(
    robust_options = list(type = "huber", max_iter = 3, scale_scope = "global"),
    ar_options = list(struct = "iid")
  )

  small <- process_run_integrated(X, Y, cfg, sigma_fixed = 0.2)
  big <- process_run_integrated(X, Y, cfg, sigma_fixed = 10)

  expect_false(isTRUE(all.equal(small$betas, big$betas, tolerance = 1e-8)))
})

test_that("solve_glm_core RSS-only path matches fitted-path results", {
  set.seed(1011)
  n <- 140
  p <- 5
  v <- 18
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rnorm(n * v), n, v)
  ctx <- glm_context(X = X, Y = Y, proj = .fast_preproject(X))

  fast <- solve_glm_core(ctx, return_fitted = FALSE)
  full <- solve_glm_core(ctx, return_fitted = TRUE)

  expect_null(fast$fitted)
  expect_equal(fast$rss, full$rss, tolerance = 1e-10)
  expect_equal(fast$sigma2, full$sigma2, tolerance = 1e-10)
  expect_equal(fast$betas, full$betas, tolerance = 1e-10)
})

test_that("fmri_lm_control accepts voxel and local robust scale scopes", {
  cfg_voxel <- fmri_lm_control(
    robust_options = list(type = "huber", scale_scope = "voxel"),
    ar_options = list(struct = "iid")
  )
  expect_equal(cfg_voxel$robust$scale_scope, "voxel")

  cfg_local <- fmri_lm_control(
    robust_options = list(type = "huber", scale_scope = "local"),
    ar_options = list(struct = "iid")
  )
  expect_equal(cfg_local$robust$scale_scope, "voxel")
})

test_that("process_run_integrated reports fixed sigma when sigma_fixed is supplied", {
  set.seed(1010)
  n <- 100
  p <- 4
  v <- 16
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rt(n * v, df = 4), n, v)

  cfg <- fmri_lm_control(
    robust_options = list(type = "huber", max_iter = 4, scale_scope = "global"),
    ar_options = list(struct = "iid")
  )

  out <- process_run_integrated(X, Y, cfg, sigma_fixed = 2.5)
  expect_equal(out$sigma, rep(2.5, v), tolerance = 1e-12)
})

test_that("voxel robust scale scope produces voxel-specific sigma", {
  set.seed(1012)
  n <- 120
  p <- 4
  v <- 12
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  scales <- seq(0.5, 3.0, length.out = v)
  E <- sweep(matrix(rt(n * v, df = 4), n, v), 2, scales, `*`)
  Y <- X %*% B + E

  cfg_voxel <- fmri_lm_control(
    robust_options = list(type = "huber", max_iter = 4, scale_scope = "voxel"),
    ar_options = list(struct = "iid")
  )
  cfg_global <- fmri_lm_control(
    robust_options = list(type = "huber", max_iter = 4, scale_scope = "global"),
    ar_options = list(struct = "iid")
  )

  out_voxel <- process_run_integrated(X, Y, cfg_voxel)
  out_global <- process_run_integrated(X, Y, cfg_global)

  expect_gt(stats::sd(out_voxel$sigma), 0)
  expect_equal(stats::sd(out_global$sigma), 0, tolerance = 1e-12)
})

test_that("process_run_integrated honors phi_fixed in AR pipeline", {
  skip_if_not_installed("fmriAR")

  set.seed(1002)
  n <- 160
  p <- 3
  v <- 12
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)

  # Add AR(1)-like noise with strong autocorrelation.
  noise <- replicate(v, as.numeric(arima.sim(model = list(ar = 0.7), n = n)))
  Y <- X %*% B + noise

  cfg <- fmri_lm_control(
    robust_options = list(type = FALSE),
    ar_options = list(struct = "ar1", iter_gls = 2, exact_first = FALSE)
  )

  phi0 <- process_run_integrated(X, Y, cfg, phi_fixed = 0)
  phi8 <- process_run_integrated(X, Y, cfg, phi_fixed = 0.8)

  expect_false(isTRUE(all.equal(phi0$betas, phi8$betas, tolerance = 1e-8)))
})

test_that("integrated solver dfres field is consumed consistently", {
  set.seed(1003)
  n <- 90
  p <- 3
  v <- 7
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rnorm(n * v), n, v)

  cfg <- fmri_lm_control(
    robust_options = list(type = FALSE),
    ar_options = list(struct = "iid")
  )

  raw <- solve_integrated_glm(X, Y, cfg)
  out <- process_run_integrated(X, Y, cfg)

  expect_true("dfres" %in% names(raw))
  expect_equal(out$dfres, raw$dfres)
})

test_that("process_run_integrated contrasts use voxelwise sigma2", {
  set.seed(1004)
  n <- 100
  p <- 3
  v <- 10
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)

  # Strongly heteroskedastic columns to ensure sigma differs by voxel.
  col_scales <- seq(0.5, 3, length.out = v)
  E <- matrix(rnorm(n * v), n, v)
  E <- sweep(E, 2, col_scales, `*`)
  Y <- X %*% B + E

  cfg <- fmri_lm_control(
    robust_options = list(type = FALSE),
    ar_options = list(struct = "iid")
  )

  w <- 1
  attr(w, "colind") <- 2L
  conlist <- list(beta2 = w)

  out <- process_run_integrated(X, Y, cfg, conlist_weights = conlist, fconlist_weights = list())
  sigma_out <- out$contrasts[[1]]$data[[1]]$sigma

  expect_length(sigma_out, v)
  expect_gt(stats::sd(sigma_out), 0)
})

test_that("coef.fmri_lm include_baseline works on pooled multi-run fits", {
  set.seed(1005)
  n_run <- 2
  run_len <- 120
  n <- n_run * run_len
  v <- 120

  onsets <- sort(sample(seq(6, run_len - 10, by = 5), 20, replace = FALSE))
  event_table <- data.frame(
    onset = c(onsets, onsets),
    cond = factor(rep(c("A", "B"), length.out = 40)),
    run = rep(1:2, each = 20)
  )

  Xsig <- cbind(1, scale(sin(seq_len(n) / 13)))
  B <- matrix(rnorm(ncol(Xsig) * v, sd = 0.35), nrow = ncol(Xsig), ncol = v)
  Y <- Xsig %*% B + matrix(rnorm(n * v, sd = 1), n, v)

  # Add sparse spikes to exercise robust/meta path.
  spike_rows <- sample((run_len + 1):n, 6)
  Y[spike_rows, ] <- Y[spike_rows, ] + matrix(rnorm(length(spike_rows) * v, sd = 8), length(spike_rows), v)

  dset <- fmridataset::matrix_dataset(
    Y,
    TR = 1,
    run_length = rep(run_len, n_run),
    event_table = event_table
  )

  sframe <- fmrihrf::sampling_frame(rep(run_len, n_run), TR = 1)
  base <- baseline_model(basis = "poly", degree = 1, sframe = sframe)

  fit <- fmri_lm(
    onset ~ hrf(cond),
    block = ~ run,
    dataset = dset,
    baseline_model = base,
    strategy = "runwise",
    use_fast_path = TRUE,
    robust = "huber",
    progress = FALSE
  )

  expect_no_error({
    b <- coef(fit, include_baseline = TRUE)
    expect_true(is.matrix(b))
    expect_equal(length(colnames(b)), ncol(b))
  })
})

test_that("AR processing routes through fmriAR in production runwise path", {
  skip_if_not_installed("fmriAR")

  set.seed(1006)
  n <- 120
  v <- 30
  event_table <- data.frame(
    onset = c(8, 24, 40, 56, 72, 88, 104),
    cond = factor(rep("A", 7)),
    run = 1
  )
  Y <- matrix(rnorm(n * v), n, v)
  dset <- fmridataset::matrix_dataset(Y, TR = 1, run_length = n, event_table = event_table)
  sframe <- fmrihrf::sampling_frame(n, TR = 1)
  base <- baseline_model(basis = "poly", degree = 1, sframe = sframe)

  .GlobalEnv$.fit_noise_calls <- 0L
  .GlobalEnv$.whiten_calls <- 0L

  trace("fit_noise",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.fit_noise_calls <- .GlobalEnv$.fit_noise_calls + 1L),
        print = FALSE)
  trace("whiten_apply",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.whiten_calls <- .GlobalEnv$.whiten_calls + 1L),
        print = FALSE)
  on.exit({
    untrace("fit_noise", where = asNamespace("fmriAR"))
    untrace("whiten_apply", where = asNamespace("fmriAR"))
    rm(.fit_noise_calls, .whiten_calls, envir = .GlobalEnv)
  }, add = TRUE)

  fmri_lm(
    onset ~ hrf(cond),
    block = ~ run,
    dataset = dset,
    baseline_model = base,
    strategy = "runwise",
    use_fast_path = TRUE,
    ar_options = list(struct = "ar1", iter_gls = 1),
    progress = FALSE
  )

  expect_gt(.GlobalEnv$.fit_noise_calls, 0L)
  expect_gt(.GlobalEnv$.whiten_calls, 0L)
})

test_that("AR processing routes through fmriAR in production chunkwise path", {
  skip_if_not_installed("fmriAR")

  set.seed(1007)
  n <- 120
  v <- 30
  event_table <- data.frame(
    onset = c(8, 24, 40, 56, 72, 88, 104),
    cond = factor(rep("A", 7)),
    run = 1
  )
  Y <- matrix(rnorm(n * v), n, v)
  dset <- fmridataset::matrix_dataset(Y, TR = 1, run_length = n, event_table = event_table)
  sframe <- fmrihrf::sampling_frame(n, TR = 1)
  base <- baseline_model(basis = "poly", degree = 1, sframe = sframe)

  .GlobalEnv$.fit_noise_calls <- 0L
  .GlobalEnv$.whiten_calls <- 0L

  trace("fit_noise",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.fit_noise_calls <- .GlobalEnv$.fit_noise_calls + 1L),
        print = FALSE)
  trace("whiten_apply",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.whiten_calls <- .GlobalEnv$.whiten_calls + 1L),
        print = FALSE)
  on.exit({
    untrace("fit_noise", where = asNamespace("fmriAR"))
    untrace("whiten_apply", where = asNamespace("fmriAR"))
    rm(.fit_noise_calls, .whiten_calls, envir = .GlobalEnv)
  }, add = TRUE)

  fmri_lm(
    onset ~ hrf(cond),
    block = ~ run,
    dataset = dset,
    baseline_model = base,
    strategy = "chunkwise",
    nchunks = 2,
    use_fast_path = TRUE,
    ar_options = list(struct = "ar1", iter_gls = 1),
    progress = FALSE
  )

  expect_gt(.GlobalEnv$.fit_noise_calls, 0L)
  expect_gt(.GlobalEnv$.whiten_calls, 0L)
})

test_that("integrated AR+robust pipeline avoids redundant whitening passes", {
  skip_if_not_installed("fmriAR")

  set.seed(1014)
  n <- 120
  p <- 4
  v <- 12
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rnorm(n * v), n, v)

  cfg <- fmri_lm_control(
    robust_options = list(type = "huber", max_iter = 3, scale_scope = "run"),
    ar_options = list(struct = "ar1", iter_gls = 2, exact_first = FALSE)
  )

  .GlobalEnv$.whiten_calls <- 0L
  trace("whiten_apply",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.whiten_calls <- .GlobalEnv$.whiten_calls + 1L),
        print = FALSE)
  on.exit({
    untrace("whiten_apply", where = asNamespace("fmriAR"))
    rm(.whiten_calls, envir = .GlobalEnv)
  }, add = TRUE)

  out <- process_run_integrated(X, Y, cfg)

  expect_true(is.list(out))
  expect_gte(.GlobalEnv$.whiten_calls, 1L)
  expect_lte(.GlobalEnv$.whiten_calls, cfg$ar$iter_gls)
})

test_that("fixed-order AR estimation avoids legacy per-voxel fallback by default", {
  skip_if_not_installed("fmriAR")

  set.seed(1008)
  resid <- replicate(40, as.numeric(arima.sim(model = list(ar = 0.6), n = 240)))
  resid <- as.matrix(resid)

  .GlobalEnv$.legacy_ar_calls <- 0L
  trace("estimate_ar_parameters",
        where = asNamespace("fmrireg"),
        tracer = quote(.GlobalEnv$.legacy_ar_calls <- .GlobalEnv$.legacy_ar_calls + 1L),
        print = FALSE)
  old_opt <- getOption("fmrireg.ar.fixed_order_legacy_fallback")
  options(fmrireg.ar.fixed_order_legacy_fallback = FALSE)
  on.exit({
    options(fmrireg.ar.fixed_order_legacy_fallback = old_opt)
    untrace("estimate_ar_parameters", where = asNamespace("fmrireg"))
    rm(.legacy_ar_calls, envir = .GlobalEnv)
  }, add = TRUE)

  phi <- fmrireg:::.estimate_phi_fixed_order(resid, order = 1L, pooling = "global")

  expect_true(is.list(phi))
  expect_equal(length(phi), 1L)
  expect_equal(length(phi[[1]]), 1L)
  expect_gt(abs(phi[[1]][1]), 0.1)
  expect_lt(abs(phi[[1]][1]), 0.99)
  expect_equal(.GlobalEnv$.legacy_ar_calls, 0L)
})

test_that("non-voxelwise AR router defaults to zero fallback without legacy estimator", {
  skip_if_not_installed("fmriAR")

  set.seed(1009)
  resid <- rnorm(200)

  local_mocked_bindings(
    .estimate_ar_via_fmriAR = function(...) list(phi = list(numeric(0))),
    .package = "fmrireg"
  )

  .GlobalEnv$.legacy_ar_calls <- 0L
  trace("estimate_ar_parameters",
        where = asNamespace("fmrireg"),
        tracer = quote(.GlobalEnv$.legacy_ar_calls <- .GlobalEnv$.legacy_ar_calls + 1L),
        print = FALSE)
  old_opt <- getOption("fmrireg.ar.nonvoxel_legacy_fallback")
  options(fmrireg.ar.nonvoxel_legacy_fallback = FALSE)
  on.exit({
    options(fmrireg.ar.nonvoxel_legacy_fallback = old_opt)
    untrace("estimate_ar_parameters", where = asNamespace("fmrireg"))
    rm(.legacy_ar_calls, envir = .GlobalEnv)
  }, add = TRUE)

  phi <- fmrireg:::.estimate_ar_parameters_routed(resid, ar_order = 2L)

  expect_equal(phi, c(0, 0))
  expect_equal(.GlobalEnv$.legacy_ar_calls, 0L)
})

test_that("non-voxelwise AR router can opt in to legacy fallback", {
  skip_if_not_installed("fmriAR")

  set.seed(1013)
  resid <- rnorm(200)

  local_mocked_bindings(
    .estimate_ar_via_fmriAR = function(...) list(phi = list(numeric(0))),
    .package = "fmrireg"
  )

  .GlobalEnv$.legacy_ar_calls <- 0L
  trace("estimate_ar_parameters",
        where = asNamespace("fmrireg"),
        tracer = quote(.GlobalEnv$.legacy_ar_calls <- .GlobalEnv$.legacy_ar_calls + 1L),
        print = FALSE)
  old_opt <- getOption("fmrireg.ar.nonvoxel_legacy_fallback")
  options(fmrireg.ar.nonvoxel_legacy_fallback = TRUE)
  on.exit({
    options(fmrireg.ar.nonvoxel_legacy_fallback = old_opt)
    untrace("estimate_ar_parameters", where = asNamespace("fmrireg"))
    rm(.legacy_ar_calls, envir = .GlobalEnv)
  }, add = TRUE)

  phi <- fmrireg:::.estimate_ar_parameters_routed(resid, ar_order = 2L)

  expect_length(phi, 2L)
  expect_gt(.GlobalEnv$.legacy_ar_calls, 0L)
})

test_that("iterative AR fmriAR path uses subset whitening before final full pass", {
  skip_if_not_installed("fmriAR")

  set.seed(1015)
  n <- 140
  p <- 4
  v <- 96
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- matrix(rnorm(p * v), nrow = p, ncol = v)
  Y <- X %*% B + matrix(rnorm(n * v), n, v)

  cfg <- list(
    struct = "ar1",
    iter_gls = 3L,
    min_iter = 99L,
    exact_first = FALSE,
    global = FALSE,
    voxelwise = FALSE,
    estimation_max_vox = 16L
  )

  .GlobalEnv$.whiten_ncol <- integer()
  trace("whiten_apply",
        where = asNamespace("fmriAR"),
        tracer = quote(.GlobalEnv$.whiten_ncol <- c(.GlobalEnv$.whiten_ncol, ncol(Y))),
        print = FALSE)
  on.exit({
    untrace("whiten_apply", where = asNamespace("fmriAR"))
    rm(.whiten_ncol, envir = .GlobalEnv)
  }, add = TRUE)

  out <- fmrireg:::.iterative_ar_gls_via_fmriAR(
    X = X,
    Y = Y,
    cfg = cfg,
    run_indices = list(seq_len(n)),
    max_iter = 3L,
    tol = 0
  )

  expect_true(is.list(out))
  expect_gte(length(.GlobalEnv$.whiten_ncol), 2L)
  expect_true(any(.GlobalEnv$.whiten_ncol < v))
  expect_equal(tail(.GlobalEnv$.whiten_ncol, 1L), v)
  expect_true(all(head(.GlobalEnv$.whiten_ncol, -1L) <= cfg$estimation_max_vox))
})

test_that("process_chunk robust-only path preserves run row mapping after in-place assembly", {
  set.seed(1016)
  n <- 10
  v <- 5
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * v), nrow = n, ncol = v)

  run_rows <- list(1:4, 5:10)
  w1 <- seq(0.7, 1.0, length.out = length(run_rows[[1]]))
  w2 <- seq(1.0, 0.6, length.out = length(run_rows[[2]]))

  precomp <- list(
    ar_modeling = FALSE,
    run_info = list(
      list(sqrtw = w1, censor = NULL),
      list(sqrtw = w2, censor = NULL)
    ),
    run_row_inds = run_rows,
    X_global = X,
    proj_global = .fast_preproject(X),
    ar_order = 0L
  )

  cfg <- fmri_lm_control(
    robust_options = list(type = "huber"),
    ar_options = list(struct = "iid")
  )

  got <- process_chunk(Y, precomp, cfg)

  Y_expected <- Y
  Y_expected[run_rows[[1]], ] <- Y[run_rows[[1]], ] * w1
  Y_expected[run_rows[[2]], ] <- Y[run_rows[[2]], ] * w2
  expected <- solve_glm_core(glm_context(X = X, Y = Y_expected, proj = precomp$proj_global))

  expect_equal(got$betas, expected$betas, tolerance = 1e-10)
  expect_equal(got$rss, expected$rss, tolerance = 1e-10)
  expect_equal(got$sigma2, expected$sigma2, tolerance = 1e-10)
})
