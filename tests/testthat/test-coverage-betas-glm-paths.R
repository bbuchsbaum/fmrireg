# High-signal coverage for estimate_betas.matrix_dataset and glm_ols/glm_lss.

make_beta_fixture <- function(n_time = 80L, n_vox = 4L, n_events = 8L, seed = 42L) {
  set.seed(seed)
  onsets <- seq(5, n_time - 10, length.out = n_events)
  etab <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = n_events)),
    nuisance = rnorm(n_events),
    run = 1L
  )
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  list(
    dataset = matrix_dataset(Y, TR = 1, run_length = n_time, event_table = etab),
    events = etab,
    Y = Y
  )
}

test_that("estimate_betas.matrix_dataset covers ols/lss/mixed and fixed effects", {
  fx <- make_beta_fixture()
  bmod <- baseline_model(basis = "constant", sframe = fx$dataset$sampling_frame)

  ols <- estimate_betas(
    fx$dataset,
    fixed = onset ~ hrf(nuisance),
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "ols",
    basemod = bmod,
    progress = FALSE
  )
  expect_s3_class(ols, "fmri_betas")
  expect_equal(ncol(ols$betas_ran), ncol(fx$Y))
  expect_true(is.matrix(ols$betas_fixed))
  expect_equal(ncol(ols$betas_fixed), ncol(fx$Y))
  expect_true(NCOL(ols$design_ran) >= 1L && NROW(ols$design_ran) >= 1L)
  expect_true(NCOL(ols$design_fixed) >= 1L && NROW(ols$design_fixed) >= 1L)
  expect_true(NCOL(ols$design_base) >= 1L && NROW(ols$design_base) >= 1L)

  lss <- estimate_betas(
    fx$dataset,
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "lss",
    progress = FALSE
  )
  expect_s3_class(lss, "fmri_betas")
  expect_equal(ncol(lss$betas_ran), ncol(fx$Y))
  expect_null(lss$betas_fixed)

  mixed <- estimate_betas(
    fx$dataset,
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "mixed",
    progress = FALSE
  )
  expect_s3_class(mixed, "fmri_betas")
  expect_equal(ncol(mixed$betas_ran), ncol(fx$Y))
})

test_that("glm_ols and glm_lss happy paths plus validation branches", {
  fx <- make_beta_fixture(n_time = 70L, n_vox = 3L, n_events = 6L, seed = 7L)
  emod <- event_model(
    onset ~ hrf(condition),
    data = fx$events,
    block = ~ run,
    sampling_frame = fx$dataset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = fx$dataset$sampling_frame)

  ols_obj <- glm_ols(
    fx$dataset, emod, fmrihrf::HRF_SPMG1,
    basemod = bmod, block = ~ run, progress = FALSE
  )
  expect_s3_class(ols_obj, "fmri_betas")
  expect_true(is.matrix(ols_obj$betas_ran))

  ols_str <- glm_ols(
    fx$dataset, emod, "HRF_SPMG1",
    block = ~ run, progress = FALSE
  )
  expect_equal(dim(ols_str$betas_ran), dim(ols_obj$betas_ran))

  expect_error(glm_ols(list(), emod, fmrihrf::HRF_SPMG1), "matrix_dataset")
  expect_error(glm_ols(fx$dataset, list(), fmrihrf::HRF_SPMG1), "event_model")
  expect_error(glm_ols(fx$dataset, emod, "not_a_basis"), "Unknown HRF")
  expect_error(glm_ols(fx$dataset, emod, 123), "HRF object")

  lss <- suppressWarnings(glm_lss(
    fx$dataset, emod, fmrihrf::HRF_SPMG1,
    block = ~ run, progress = FALSE
  ))
  expect_s3_class(lss, "fmri_betas")
  expect_true(is.matrix(lss$betas_ran))

  expect_warning(
    glm_lss(
      fx$dataset, emod, "HRF_SPMG1",
      block = ~ run, use_cpp = TRUE, progress = FALSE
    ),
    "C\\+\\+-optimized LSS|retired"
  )

  expect_error(glm_lss(list(), emod, fmrihrf::HRF_SPMG1), "matrix_dataset")
  expect_error(glm_lss(fx$dataset, list(), fmrihrf::HRF_SPMG1), "event_model")
  expect_error(glm_lss(fx$dataset, emod, "bad_basis"), "Unknown HRF")
  expect_error(glm_lss(fx$dataset, emod, TRUE), "HRF object")
})

test_that("mixed_betas validation and inject_basis rewrite formulas", {
  X <- cbind(1, rnorm(24), rnorm(24))
  y <- rnorm(24)

  expect_error(fmrireg:::mixed_betas(matrix(0, 0, 0), 1, 1L, NULL), "zero rows")
  expect_error(fmrireg:::mixed_betas(X, y, integer(0), NULL), "No random")
  expect_error(fmrireg:::mixed_betas(X, y, 1:2, 99L), "Index out of bounds")
  expect_error(fmrireg:::mixed_betas(X, y[1:10], 1:2, NULL), "same number of rows")
  expect_error(fmrireg:::mixed_betas(X, y, 1:2, NULL, solver = "nope"), "NULL or a function")
  expect_error(
    fmrireg:::mixed_betas(X, y, c(1.5, 2), NULL),
    "positive integers"
  )

  ok <- fmrireg:::mixed_betas(X, y, ran_ind = 2:3, fixed_ind = 1L)
  expect_equal(length(ok), 3L)

  ok_nofixed <- fmrireg:::mixed_betas(X, y, ran_ind = 2:3, fixed_ind = NULL)
  expect_equal(length(ok_nofixed), 2L)

  form <- onset ~ hrf(condition) + trialwise(amp)
  rewritten <- fmrireg:::inject_basis(form, fmrihrf::HRF_SPMG1)
  expect_true(inherits(rewritten, "formula"))
  rhs <- paste(deparse(rewritten[[3]]), collapse = " ")
  expect_true(grepl("basis", rhs))
})
