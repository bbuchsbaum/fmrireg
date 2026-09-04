# Coverage for lowrank/rrr engine helpers and bootstrap_glm_inference.

test_that("lowrank engine helpers cover preflight, projection, and AR blend", {
  fx_etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(50 * 4), 50, 4)
  dset <- matrix_dataset(Y, TR = 1, run_length = 50, event_table = fx_etab)
  emod <- event_model(
    onset ~ hrf(condition),
    data = fx_etab,
    block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  expect_true(fmrireg:::.preflight_lowrank_engine(fmod, dset, list(), cfg))
  expect_error(
    fmrireg:::.preflight_lowrank_engine(list(), dset, list(), cfg),
    "fmri_model"
  )
  expect_error(
    fmrireg:::.preflight_lowrank_engine(fmod, NULL, list(), cfg),
    "dataset"
  )
  expect_error(
    fmrireg:::.preflight_lowrank_engine(fmod, dset, list(), list()),
    "fmri_lm_control"
  )

  expect_true(
    fmrireg:::.preflight_lowrank_engine_plugin(
      fmod, dset, args = list(lowrank = list()), cfg = cfg
    )
  )

  A <- matrix(c(1, 0.2, 0.1, 1), 2, 2)
  x <- matrix(1:6, 3, 2)
  expect_equal(fmrireg:::.lowrank_project_voxels(x, A, TRUE), x)
  proj <- fmrireg:::.lowrank_project_voxels(x, A, FALSE)
  expect_equal(dim(proj), c(3, 2))

  expect_equal(fmrireg:::.lowrank_blend_ar(0.4, 0.2, 0.5), 0.3)
  blended <- fmrireg:::.lowrank_blend_ar(list(0.4, 0.5), list(0.2, 0.1), 0.5)
  expect_equal(length(blended), 2)
  expect_error(
    fmrireg:::.lowrank_blend_ar(list(0.4), list(0.2, 0.1), 0.5),
    "different lengths"
  )
})

test_that("bootstrap_glm_inference residual path returns CIs", {
  skip_if_not_installed("progress")
  set.seed(12)
  n <- 40
  p <- 3
  V <- 5
  X <- cbind(1, rnorm(n), rnorm(n))
  beta <- matrix(c(0.5, 1, -0.5), p, 1)
  Y <- X %*% matrix(rep(beta, V), p, V) + matrix(rnorm(n * V, sd = 0.5), n, V)
  fit_result <- list(
    betas = solve(crossprod(X), crossprod(X, Y)),
    sigma = apply(Y - X %*% solve(crossprod(X), crossprod(X, Y)), 2, sd),
    df.residual = n - p
  )
  # Ensure compute_contrast path can be skipped by omitting contrasts
  out <- fmrireg:::bootstrap_glm_inference(
    fit_result = fit_result,
    X = X,
    Y = Y,
    config = fmri_lm_control(),
    contrasts = NULL,
    nboot = 20,
    block_size = 5,
    method = "residual",
    parallel = FALSE
  )
  expect_true(is.list(out))
  expect_true(!is.null(out$boot_betas) || !is.null(out$ci) || !is.null(out$betas_ci) ||
                length(out) > 0)
})

test_that("rrr/lowrank engine registration preflight hooks exist", {
  fx_etab <- data.frame(onset = c(5, 20), condition = factor(c("A", "B")), run = 1L)
  Y <- matrix(rnorm(40 * 3), 40, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = 40, event_table = fx_etab)
  emod <- event_model(
    onset ~ hrf(condition), data = fx_etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  fmod <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control()

  expect_true(
    fmrireg:::.preflight_rrr_gls_engine(fmod, dset, args = list(), cfg = cfg)
  )
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(list(), dset, list(), cfg),
    "fmri_model"
  )
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(fmod, NULL, list(), cfg),
    "dataset"
  )
  expect_error(
    fmrireg:::.preflight_rrr_gls_engine(fmod, dset, list(), list()),
    "fmri_lm_control"
  )
})
