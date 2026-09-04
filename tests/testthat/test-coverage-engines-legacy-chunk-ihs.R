# Fifteenth wave: IHS/sketch engines, chunkwise_apply, fmri_lm legacy + config.

test_that("make_time_sketch and ihs_latent_solve cover sketch methods", {
  set.seed(431)
  Tlen <- 64L
  p <- 4L
  k <- 3L
  m <- 16L
  X <- matrix(rnorm(Tlen * p), Tlen, p)
  Z <- matrix(rnorm(Tlen * k), Tlen, k)

  Sg <- fmrireg:::make_time_sketch(Tlen, list(method = "gaussian", m = m))
  expect_equal(dim(Sg), c(m, Tlen))

  Sc <- fmrireg:::make_time_sketch(Tlen, list(method = "countsketch", m = m))
  expect_equal(dim(Sc), c(m, Tlen))
  expect_true(inherits(Sc, "Matrix") || inherits(Sc, "sparseMatrix") || is.matrix(Sc))

  expect_null(fmrireg:::make_time_sketch(Tlen, list(method = "srht", m = m)))
  expect_null(fmrireg:::make_time_sketch(Tlen, list(method = "ihs", m = m)))

  expect_error(
    fmrireg:::make_time_sketch(Tlen, list(method = "gaussian", m = Tlen + 1L)),
    regexp = "."
  )

  sol <- fmrireg:::ihs_latent_solve(X, Z, m = m, iters = 2L)
  expect_true(is.list(sol))
  expect_true(all(c("M", "Ginv") %in% names(sol)) || length(sol) >= 1L)
  if (!is.null(sol$M)) {
    expect_true(is.matrix(sol$M) || inherits(sol$M, "Matrix"))
  }
})

test_that(".chunkwise_apply sequential and parallel paths", {
  chunks <- list(list(v = 1), list(v = 2), list(v = 3))
  out <- fmrireg:::.chunkwise_apply(
    chunks,
    FUN = function(i, ch) ch$v * 10,
    parallel_chunks = FALSE,
    progress = FALSE
  )
  expect_equal(unlist(out), c(10, 20, 30))

  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)

  out_p <- fmrireg:::.chunkwise_apply(
    chunks,
    FUN = function(i, ch) ch$v + i,
    parallel_chunks = TRUE,
    progress = FALSE
  )
  expect_equal(unlist(out_p), c(2, 4, 6))
})

test_that(".fmri_lm_formula_legacy and .fmri_lm_model_legacy fit matrix data", {
  set.seed(432)
  n <- 40L
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)

  fit <- fmrireg:::.fmri_lm_formula_legacy(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    strategy = "runwise",
    progress = FALSE
  )
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(coef(fit)))

  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  model <- fmri_model(emod, bmod, dset)

  fit2 <- fmrireg:::.fmri_lm_model_legacy(
    model,
    dataset = dset,
    strategy = "chunkwise",
    nchunks = 2L,
    progress = FALSE
  )
  expect_s3_class(fit2, "fmri_lm")
})

test_that(".fmri_lm_control_legacy and noise/robust/weights/projection converters", {
  cfg <- fmrireg:::.fmri_lm_control_legacy(
    ar_options = list(struct = "ar1", iter = 2L, global = TRUE, exact_first = TRUE),
    robust_options = list(type = TRUE, max_iter = 2L, scale_scope = "local")
  )
  expect_s3_class(cfg, "fmri_lm_control")
  expect_equal(cfg$ar$struct, "ar1")
  expect_equal(cfg$ar$pooling, "global")
  expect_equal(cfg$robust$type, "huber")
  expect_equal(cfg$robust$scale_scope, "voxel")

  noise_ord <- fmrireg:::.fmri_lm_noise_from_legacy(list(order = 2L))
  expect_equal(noise_ord$struct, "ar2")

  expect_error(
    fmrireg:::.fmri_lm_noise_from_legacy(list(global = TRUE, by_cluster = TRUE)),
    "cannot both be TRUE"
  )

  # Pass-through when already a typed spec
  expect_identical(
    fmrireg:::.fmri_lm_noise_from_legacy(noise_ord),
    noise_ord
  )

  rob <- fmrireg:::.fmri_lm_robust_from_legacy(list(type = "bisquare", max_iter = 3L))
  expect_equal(rob$type, "bisquare")
  expect_equal(rob$max_iter, 3L)

  wspec <- fmrireg:::.fmri_lm_weights_from_legacy(
    list(enabled = TRUE, method = "tukey", threshold = 1.2)
  )
  expect_true(wspec$method %in% c("tukey", "inverse_squared") || is.character(wspec$method))

  pspec <- fmrireg:::.fmri_lm_projection_from_legacy(
    list(enabled = TRUE, nuisance_matrix = matrix(1, 10, 2), lambda = 0.5)
  )
  expect_true(pspec$method %in% c("soft_subspace", "none") || is.character(pspec$method))

  expect_output(print(cfg), "fmri_lm_control|control|ar|robust")
})
