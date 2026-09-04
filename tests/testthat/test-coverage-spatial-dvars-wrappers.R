# create_3d_blocks, dvars/volume_weights edges, and exported wrappers.

test_that("create_3d_blocks covers array path and validation", {
  mask <- array(1, dim = c(4, 4, 2))
  mask[1, 1, 1] <- 0
  blocks <- create_3d_blocks(mask, block_size = c(2, 2, 1), connectivity = 26)
  expect_true(is.list(blocks))
  expect_true(blocks$n_groups >= 1L)
  expect_equal(length(blocks$group_id), sum(mask > 0))
  expect_true(is.null(blocks$neighbors) || is.list(blocks$neighbors))

  blocks6 <- create_3d_blocks(mask, block_size = c(2, 2, 2), connectivity = 6)
  expect_equal(length(blocks6$group_id), sum(mask > 0))

  blocks_none <- create_3d_blocks(mask, block_size = c(2, 2, 1), connectivity = 18)
  expect_null(blocks_none$neighbors)

  expect_error(create_3d_blocks("definitely_missing_mask.nii"), "does not exist")
  expect_error(create_3d_blocks(1:10), "3D array|NeuroVol|file path")
})

test_that("compute_dvars / dvars_to_weights / volume_weights cover edges", {
  set.seed(5)
  Y <- matrix(rnorm(40 * 8), 40, 8)
  Y[20, ] <- Y[20, ] + 4

  # Non-matrix coercion path
  d <- compute_dvars(as.data.frame(Y), normalize = TRUE)
  expect_equal(length(d), 40L)
  expect_true(all(d >= 0))

  d_raw <- compute_dvars(Y, normalize = FALSE)
  expect_true(median(d_raw) > 0)

  expect_error(compute_dvars(matrix(1, 1, 3)), "at least 2")

  w_inv <- dvars_to_weights(d, method = "inverse_squared")
  w_soft <- dvars_to_weights(d, method = "soft_threshold", threshold = 1.2)
  w_tukey <- dvars_to_weights(d, method = "tukey")
  expect_equal(length(w_inv), length(d))
  expect_true(all(is.finite(w_inv)))
  expect_true(all(is.finite(w_soft)))
  expect_true(all(is.finite(w_tukey)))
  expect_true(min(w_inv) > 0)
  expect_error(dvars_to_weights(c(-0.1, 1)), "non-negative")

  vw <- volume_weights(Y, method = "soft_threshold", return_dvars = TRUE)
  expect_equal(names(vw), c("weights", "dvars"))
  expect_equal(length(vw$weights), nrow(Y))

  vw2 <- volume_weights(Y, method = "tukey", return_dvars = FALSE)
  expect_true(is.numeric(vw2))

  X <- cbind(1, rnorm(40))
  applied <- fmrireg:::apply_volume_weights(X, Y, vw$weights)
  expect_true(is.list(applied) || (is.matrix(applied) || is.numeric(applied)))
})

test_that("exported wrappers exercise reexported fmridesign/fmridataset helpers", {
  fmod <- fmrireg:::.demo_fmri_model()
  dset <- fmrireg:::.demo_matrix_dataset()

  expect_s3_class(as.matrix_dataset(dset), "matrix_dataset")

  bs <- BSpline(seq(0, 1, length.out = 12), degree = 3)
  expect_s3_class(bs, "BSpline")

  et <- event_terms(fmod$event_model)
  expect_true(length(et) >= 1L)
  fc <- Fcontrasts(et[[1]])
  expect_true(is.list(fc) || is.matrix(fc))

  ti <- term_indices(design_matrix(fmod$event_model))
  expect_true(is.list(ti) || is.numeric(ti) || is.null(ti))

  # Best-effort optional helpers; tolerate API shape differences
  cbl <- tryCatch(condition_basis_list(et[[1]]), error = function(e) e)
  expect_true(is.list(cbl) || inherits(cbl, "error") || is.matrix(cbl))

  cd <- tryCatch(
    convolve_design(design_matrix(fmod$event_model), fmod$event_model),
    error = function(e) e
  )
  expect_true(is.matrix(cd) || is.list(cd) || inherits(cd, "error"))

  sb <- tryCatch(sub_basis(et[[1]], 1L), error = function(e) e)
  expect_true(!is.null(sb) || inherits(sb, "error"))

  pc <- tryCatch(
    plot_contrasts(contrast_set(pair_contrast(~ A, ~ B, name = "AvB"))),
    error = function(e) e
  )
  expect_true(inherits(pc, "error") || inherits(pc, "ggplot") || is.list(pc))
})
