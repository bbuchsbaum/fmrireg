# fmrilm.R: remaining robust/AR/preprocessing normalization and config attach.

test_that(".fmri_lm_normalize_robust_options resolves aliases and conflicts", {
  expect_equal(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = TRUE)$type,
    "huber"
  )
  expect_equal(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = FALSE)$type,
    FALSE
  )
  expect_equal(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = "bisquare")$type,
    "bisquare"
  )
  expect_equal(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = "false")$type,
    FALSE
  )
  expect_equal(
    fmrireg:::.fmri_lm_normalize_robust_options(
      robust = NULL, robust_psi = "huber", robust_max_iter = 5L,
      robust_scale_scope = "global"
    )$max_iter,
    5L
  )

  expect_error(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = c(TRUE, FALSE)),
    "single non-missing"
  )
  expect_error(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = 3),
    "TRUE, FALSE"
  )
  expect_error(
    fmrireg:::.fmri_lm_normalize_robust_options(robust = "tukey"),
    "TRUE, FALSE"
  )
  expect_error(
    fmrireg:::.fmri_lm_normalize_robust_options(
      robust = TRUE, robust_options = list(type = "bisquare")
    ),
    "Conflicting robust"
  )

  # engine_robust_options merges
  merged <- fmrireg:::.fmri_lm_normalize_robust_options(
    robust = NULL,
    robust_options = list(),
    engine_robust_options = list(type = "huber", max_iter = 3L)
  )
  expect_equal(merged$type, "huber")
  expect_equal(merged$max_iter, 3L)
})

test_that(".fmri_lm_normalize_ar_options maps order and rejects conflicts", {
  ar <- fmrireg:::.fmri_lm_normalize_ar_options(
    ar_options = list(order = 0L)
  )
  expect_equal(ar$struct, "iid")

  ar1 <- fmrireg:::.fmri_lm_normalize_ar_options(ar_options = list(order = 1L))
  expect_equal(ar1$struct, "ar1")
  ar2 <- fmrireg:::.fmri_lm_normalize_ar_options(ar_options = list(order = 2L))
  expect_equal(ar2$struct, "ar2")
  arp <- fmrireg:::.fmri_lm_normalize_ar_options(ar_options = list(order = 5L))
  expect_equal(arp$struct, "arp")
  expect_equal(arp$p, 5L)

  # Shorthand aliases
  aliased <- fmrireg:::.fmri_lm_normalize_ar_options(
    cor_struct = "ar1", cor_iter = 2L, cor_global = TRUE,
    ar1_exact_first = TRUE, ar_p = NULL, ar_voxelwise = FALSE
  )
  expect_equal(aliased$struct, "ar1")
  expect_equal(aliased$iter_gls, 2L)
  expect_true(isTRUE(aliased$global))
  expect_true(isTRUE(aliased$exact_first))
  expect_false(isTRUE(aliased$voxelwise))

  expect_error(
    fmrireg:::.fmri_lm_normalize_ar_options(
      ar_options = list(struct = "ar2"), cor_struct = "ar1"
    ),
    "Conflicting AR"
  )
})

test_that(".fmri_lm_normalize_preprocessing_options and build/attach config", {
  prep <- fmrireg:::.fmri_lm_normalize_preprocessing_options(
    volume_weights = "dvars",
    soft_subspace_options = NULL,
    nuisance_projection = matrix(rnorm(20), 10, 2)
  )
  expect_true(isTRUE(prep$volume_weights_options$enabled))
  expect_equal(prep$volume_weights_options$method, "dvars")
  expect_true(isTRUE(prep$soft_subspace_options$enabled))
  expect_true(is.matrix(prep$soft_subspace_options$nuisance_matrix))

  prep2 <- fmrireg:::.fmri_lm_normalize_preprocessing_options(
    nuisance_projection = "mask.nii"
  )
  expect_equal(prep2$soft_subspace_options$nuisance_mask, "mask.nii")

  cfg <- fmrireg:::.fmri_lm_build_config(
    robust = "huber",
    cor_struct = "ar1",
    ar_voxelwise = FALSE
  )
  expect_s3_class(cfg, "fmri_lm_control")
  expect_true(!is.null(cfg$robust$type) || !is.null(cfg$ar$struct))

  fit <- structure(list(ok = TRUE), class = "fmri_lm")
  attached <- fmrireg:::.fmri_lm_attach_config_metadata(fit, cfg, cfg)
  expect_s3_class(attr(attached, "requested_config"), "fmri_lm_control")
  expect_s3_class(attr(attached, "executed_config"), "fmri_lm_control")
  expect_identical(attr(attached, "config"), cfg)
})

test_that(".fmri_lm_resolve_execution canonical vs legacy and conflicts", {
  cfg <- fmri_lm_control()
  legacy <- fmrireg:::.fmri_lm_resolve_execution(
    cfg, compute = NULL,
    strategy = "chunkwise", nchunks = 2L, use_fast_path = TRUE,
    progress = FALSE, parallel_voxels = FALSE, parallel_chunks = FALSE
  )
  expect_equal(legacy$strategy, "chunkwise")
  expect_equal(legacy$nchunks, 2L)
  expect_true(isTRUE(legacy$use_fast_path))

  compute <- compute_spec(
    voxel_chunks = 3L, backend = "matrix", parallel = "chunks", progress = FALSE
  )
  canon <- fmrireg:::.fmri_lm_resolve_execution(
    cfg, compute = compute,
    strategy = "runwise", nchunks = 1L, use_fast_path = FALSE,
    progress = TRUE, parallel_voxels = FALSE, parallel_chunks = FALSE
  )
  expect_equal(canon$nchunks, 3L)
  expect_true(isTRUE(canon$use_fast_path))
  expect_true(isTRUE(canon$parallel_chunks))

  expect_error(
    fmrireg:::.fmri_lm_resolve_execution(
      cfg, compute = list(fake = TRUE),
      strategy = "chunkwise", nchunks = 1L, use_fast_path = TRUE,
      progress = FALSE, parallel_voxels = FALSE, parallel_chunks = FALSE
    ),
    "compute_spec"
  )

  expect_error(
    fmrireg:::.fmri_lm_resolve_execution(
      cfg, compute = compute,
      strategy = "chunkwise", nchunks = 1L, use_fast_path = TRUE,
      progress = FALSE, parallel_voxels = FALSE, parallel_chunks = FALSE,
      legacy_supplied = c(strategy = TRUE)
    ),
    "cannot be combined"
  )
})
