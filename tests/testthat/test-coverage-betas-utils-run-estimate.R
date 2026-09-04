# Fifteenth wave: beta_utils + run_estimate_betas OLS/mixed (no progress bar).

make_beta_fixture <- function(n = 48L, V = 3L, seed = 411L) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  bmod <- baseline_model("constant", sframe = dset$sampling_frame)
  bdes <- fmrireg:::gen_beta_design(
    fixed = NULL,
    ran = onset ~ hrf(condition),
    block = ~ run,
    bmod = bmod,
    dset = dset
  )
  list(dset = dset, bdes = bdes, Y = Y)
}

test_that("build_design_data and map_voxels cover fixed/null and progress=FALSE", {
  fx <- make_beta_fixture()
  xd <- fmrireg:::build_design_data(fx$bdes)
  expect_true(is.matrix(xd$X) || inherits(xd$X, "Matrix"))
  expect_true(is.matrix(xd$Base))
  expect_equal(nrow(xd$X), nrow(xd$Base))
  expect_true(all(is.finite(xd$X)))

  # With fixed effects, X binds ran + fixed
  bdes_fix <- fx$bdes
  bdes_fix$dmat_fixed <- matrix(rnorm(nrow(fx$Y)), nrow(fx$Y), 1)
  colnames(bdes_fix$dmat_fixed) <- "fixed1"
  xd2 <- fmrireg:::build_design_data(bdes_fix)
  expect_equal(ncol(xd2$X), ncol(as.matrix(fx$bdes$dmat_ran)) + 1L)

  # NA coercion path
  bdes_na <- fx$bdes
  dmat <- as.matrix(bdes_na$dmat_ran)
  dmat[1, 1] <- NA_real_
  bdes_na$dmat_ran <- dmat
  bdes_na$dmat_fixed <- NULL
  xd3 <- fmrireg:::build_design_data(bdes_na)
  expect_false(anyNA(xd3$X))

  vecs <- fmrireg:::masked_vectors(fx$dset)
  expect_true(length(vecs) >= 1L || inherits(vecs, "NeuroVecSeq") || is.list(vecs) ||
                inherits(vecs, "vectorseq") || length(as.list(vecs)) >= 1L)

  mapped <- fmrireg:::map_voxels(vecs, function(v) mean(v), .progress = FALSE)
  expect_true(is.matrix(mapped) || is.numeric(mapped))
  expect_equal(NCOL(mapped), ncol(fx$Y))
})

test_that("run_estimate_betas OLS and mixed return beta matrices", {
  fx <- make_beta_fixture(seed = 412L)

  ols <- fmrireg:::run_estimate_betas(
    fx$bdes, fx$dset, method = "ols", block = ~ run, progress = FALSE
  )
  expect_true(is.list(ols))
  expect_true(is.matrix(ols$beta_matrix))
  expect_equal(ncol(ols$beta_matrix), ncol(fx$Y))
  expect_equal(nrow(ols$beta_matrix), length(fx$bdes$ran_ind))
  expect_null(ols$estimated_hrf)

  mixed <- fmrireg:::run_estimate_betas(
    fx$bdes, fx$dset, method = "mixed", block = ~ run, progress = FALSE
  )
  expect_true(is.matrix(mixed$beta_matrix))
  expect_equal(ncol(mixed$beta_matrix), ncol(fx$Y))
  expect_true(all(is.finite(mixed$beta_matrix) | is.na(mixed$beta_matrix)))

  expect_error(
    fmrireg:::run_estimate_betas(
      fx$bdes, fx$dset, method = "not_a_method", block = ~ run, progress = FALSE
    ),
    "mixed|lss|ols|arg"
  )
})

test_that("estimate_betas.matrix_dataset OLS with progress=FALSE", {
  fx <- make_beta_fixture(seed = 413L)
  out <- estimate_betas(
    fx$dset,
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "ols",
    progress = FALSE
  )
  expect_s3_class(out, "fmri_betas")
  expect_true(is.matrix(out$betas_ran))
  expect_equal(ncol(out$betas_ran), ncol(fx$Y))
  expect_true(!is.null(out$design_ran))
  expect_true(!is.null(out$design_base))
})
