test_that("typed specs expose strict schemas", {
  expect_s3_class(estimation_spec(), "fmri_lm_estimation_spec")
  expect_s3_class(noise_spec(), "fmri_lm_noise_spec")
  expect_s3_class(robust_spec(), "fmri_lm_robust_spec")
  expect_s3_class(variance_spec(), "fmri_lm_variance_spec")
  expect_s3_class(weights_spec(), "fmri_lm_weights_spec")
  expect_s3_class(projection_spec(), "fmri_lm_projection_spec")
  expect_s3_class(compute_spec(), "fmri_lm_compute_spec")

  expect_error(noise_spec(typo_key = 1), "unused argument")
  expect_error(robust_spec(type = c("none", "huber")), "length 1")
  expect_error(variance_spec(max_lag = -1), "must be >= 0")
  expect_error(compute_spec(voxel_chunks = 1.5), "must be an integer")
  expect_error(projection_spec(method = "soft_subspace"), "requires")
  expect_error(noise_spec(struct = "arp"), "must be supplied")
  expect_identical(noise_spec(q = 1L)$q, 1L)
})

test_that("control rejects unknown legacy keys and canonical conflicts", {
  expect_error(fmri_lm_control(ar_options = list(typo_key = 1)),
               "Unknown ar_options key")
  expect_error(fmri_lm_control(robust_options = list(typo_key = 1)),
               "Unknown robust_options key")
  expect_error(
    fmri_lm_control(noise = noise_spec(), ar_options = list(struct = "ar1")),
    "cannot be combined"
  )
  expect_error(
    fmri_lm_control(robust = robust_spec(), robust_options = list(type = "huber")),
    "cannot be combined"
  )
})

test_that("control has one canonical vocabulary and serializes losslessly", {
  control <- fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("arp", p = 3L, pooling = "global"),
    robust = robust_spec("bisquare", max_iter = 4L),
    variance = variance_spec("model", df = "satterthwaite"),
    weights = weights_spec("tukey"),
    projection = projection_spec(),
    na_action = "propagate"
  )

  expect_s3_class(control, "fmri_lm_control")
  expect_false(inherits(control, "fmri_lm_config"))
  expect_identical(control$noise, control$ar)
  expect_identical(control$weights, control$volume_weights)
  expect_identical(control$projection, control$soft_subspace)

  encoded <- paste(capture.output(dput(control)), collapse = "\n")
  decoded <- eval(parse(text = encoded), envir = baseenv())
  expect_identical(decoded, control)
})

test_that("public fit methods expose only the redesigned boundary", {
  formula_args <- names(formals(fmrireg:::fmri_lm.formula))
  model_args <- names(formals(fmrireg:::fmri_lm.fmri_model))
  canonical <- c("control", "compute", "engine", "engine_args", "lowrank", "...")
  retired <- c(
    "strategy", "nchunks", "use_fast_path", "progress",
    "parallel_voxels", "parallel_chunks", "robust", "robust_options",
    "ar_options", "ar_voxelwise", "cor_struct", "cor_iter",
    "volume_weights", "nuisance_projection", "cfg"
  )

  expect_true(all(canonical %in% formula_args))
  expect_true(all(canonical %in% model_args))
  expect_length(intersect(formula_args, retired), 0L)
  expect_length(intersect(model_args, retired), 0L)

  control_args <- names(formals(fmri_lm_control))
  expect_true(all(c("estimation", "noise", "robust", "variance", "weights",
                    "projection", "na_action", "...") %in% control_args))
  expect_length(intersect(control_args, c("robust_options", "ar_options",
                                          "volume_weights_options",
                                          "soft_subspace_options")), 0L)
})

test_that("compute is mechanical and estimation scope chooses the estimator", {
  joint <- fmri_lm_control(estimation = estimation_spec("joint"))
  runwise <- fmri_lm_control(estimation = estimation_spec("runwise_meta"))

  joint_exec <- fmrireg:::.fmri_lm_resolve_execution(
    joint, compute_spec(voxel_chunks = 3L, backend = "reference", parallel = "chunks"),
    strategy = "runwise", nchunks = 10L, use_fast_path = TRUE,
    progress = FALSE, parallel_voxels = FALSE, parallel_chunks = FALSE
  )
  run_exec <- fmrireg:::.fmri_lm_resolve_execution(
    runwise, compute_spec(voxel_chunks = 2L, parallel = "voxels"),
    strategy = "chunkwise", nchunks = 10L, use_fast_path = TRUE,
    progress = FALSE, parallel_voxels = FALSE, parallel_chunks = FALSE
  )

  expect_identical(joint_exec$strategy, "chunkwise")
  expect_identical(joint_exec$nchunks, 3L)
  expect_false(joint_exec$use_fast_path)
  expect_true(joint_exec$parallel_chunks)
  expect_identical(run_exec$strategy, "runwise")
  expect_true(run_exec$parallel_voxels)

  expect_error(
    fmrireg:::.fmri_lm_resolve_execution(
      joint, compute_spec(), strategy = "runwise", nchunks = 10L,
      use_fast_path = TRUE, progress = FALSE, parallel_voxels = FALSE,
      parallel_chunks = FALSE, legacy_supplied = c(strategy = TRUE)
    ),
    "cannot be combined"
  )
})

test_that("canonical and translated legacy calls agree for each fit scope", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(onsets ~ hrf(condition), ~run, dataset = dset)

  canonical_joint <- suppressWarnings(fmri_lm(
    model, dataset = dset,
    control = fmri_lm_control(estimation = estimation_spec("joint")),
    compute = compute_spec(voxel_chunks = 1L)
  ))
  legacy_joint <- suppressWarnings(fmri_lm(
    model, dataset = dset, strategy = "chunkwise", nchunks = 1L
  ))

  canonical_runwise <- suppressWarnings(fmri_lm(
    model, dataset = dset,
    control = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
    compute = compute_spec(voxel_chunks = 1L)
  ))
  legacy_runwise <- suppressWarnings(fmri_lm(
    model, dataset = dset, strategy = "runwise"
  ))

  expect_equal(coef(canonical_joint), coef(legacy_joint), tolerance = 1e-10)
  expect_equal(coef(canonical_runwise), coef(legacy_runwise), tolerance = 1e-10)
  expect_identical(attr(canonical_joint, "strategy"), "chunkwise")
  expect_identical(attr(canonical_runwise, "strategy"), "runwise")
})
