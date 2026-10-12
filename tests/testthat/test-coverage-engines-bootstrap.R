# Coverage for lowrank/rrr engine helpers. (The bootstrap_glm_inference test
# was removed with that function in review pass 3, R3-5.)

test_that("lowrank engine helpers cover preflight, projection, and AR blend", {
  fx_etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(50 * 4), 50, 4)
  dset <- matrix_frame(Y, TR = 1, run_length = 50, event_table = fx_etab)
  emod <- event_model(
    onset ~ hrf(condition),
    data = fx_etab,
    block = ~ run,
    sampling_frame = frame_sframe(dset)
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = frame_sframe(dset))
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

test_that("rrr/lowrank engine registration preflight hooks exist", {
  fx_etab <- data.frame(onset = c(5, 20), condition = factor(c("A", "B")), run = 1L)
  Y <- matrix(rnorm(40 * 3), 40, 3)
  dset <- matrix_frame(Y, TR = 1, run_length = 40, event_table = fx_etab)
  emod <- event_model(
    onset ~ hrf(condition), data = fx_etab, block = ~ run,
    sampling_frame = frame_sframe(dset)
  )
  bmod <- baseline_model(basis = "constant", sframe = frame_sframe(dset))
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
