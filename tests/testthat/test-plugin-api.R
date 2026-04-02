test_that("register_basis allows plugins to supply HRF constructors", {
  basis_name <- "plugin_basis_test"
  calls <- list()
  register_basis(basis_name, function(alpha = NULL, ...) {
    calls[[length(calls) + 1L]] <<- list(alpha = alpha)
    fmrihrf::HRF_SPMG1
  })
  on.exit(rm(list = basis_name, envir = fmrireg:::.fmrireg_basis_registry), add = TRUE)

  mdl <- create_fmri_model(
    formula = onsets ~ hrf(condition, basis = basis_name, alpha = 0.42),
    block = ~run,
    dataset = .demo_matrix_dataset()
  )

  expect_s3_class(mdl, "fmri_model")
  expect_equal(length(calls), 1)
  expect_equal(calls[[1]]$alpha, 0.42)
})

test_that("engine spec accessor exposes normalized metadata", {
  latent_spec <- engine_spec("latent_sketch")
  expect_s3_class(latent_spec, "fmrireg_engine_spec")
  expect_identical(latent_spec$name, "latent_sketch")
  expect_identical(latent_spec$source, "builtin")
  expect_true("sketch" %in% latent_spec$aliases)
  expect_identical(latent_spec$strategy, "sketch")
  expect_false(latent_spec$capabilities$robust)
  expect_null(latent_spec$fit)

  specs <- engine_specs()
  expect_true(all(c("latent_sketch", "rrr_gls") %in% names(specs)))
  expect_s3_class(specs[["rrr_gls"]], "fmrireg_engine_spec")
})

test_that("engine spec print method is compact and informative", {
  out <- utils::capture.output(print(engine_spec("rrr_gls")))
  txt <- paste(out, collapse = "\n")

  expect_match(txt, "<fmrireg_engine_spec>", fixed = TRUE)
  expect_match(txt, "name: rrr_gls", fixed = TRUE)
  expect_match(txt, "source: builtin | strategy: engine", fixed = TRUE)
  expect_match(txt, "requires: event regressors", fixed = TRUE)
})

test_that("register_engine integrates with fmri_lm via fit_glm_on_transformed_series", {
  engine_name <- "plugin_engine_test"
  preflight_called <- 0L
  captured <- NULL

  register_engine(
    engine_name,
    preflight = function(model, dataset, args, cfg) {
      expect_s3_class(model, "fmri_model")
      expect_true(inherits(dataset, "matrix_dataset"))
      expect_true(is.list(args))
      expect_s3_class(cfg, "fmri_lm_config")
      preflight_called <<- preflight_called + 1L
    },
    fit = function(model, dataset, args, cfg) {
      captured <<- list(model = model, dataset = dataset, args = args)
      Y <- as.matrix(fmridataset::get_data_matrix(dataset))
      suppressWarnings(
        fit_glm_on_transformed_series(
          model,
          Y,
          cfg = cfg,
          dataset = dataset,
          engine = engine_name,
          strategy = "engine"
        )
      )
    }
  )
  on.exit(rm(list = engine_name, envir = fmrireg:::.fmrireg_engine_registry), add = TRUE)

  dset <- .demo_matrix_dataset()
  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = engine_name,
    engine_args = list(k = 3),
    plugin_engine_test = list(flag = TRUE)
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(preflight_called, 1L)
  expect_true(is.list(captured))
  expect_equal(captured$args$k, 3)
  expect_true(captured$args$flag)
  expect_identical(captured$dataset, dset)
  expect_equal(attr(fit, "engine"), engine_name)
  expect_equal(attr(fit, "strategy"), "engine")
  expect_s3_class(attr(fit, "requested_config"), "fmri_lm_config")
  expect_s3_class(attr(fit, "executed_config"), "fmri_lm_config")
  expect_equal(attr(fit, "config"), attr(fit, "executed_config"))
})

test_that("sketch engine alias dispatches to latent sketch path", {
  dset <- .demo_matrix_dataset()

  fit <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = "sketch"
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(attr(fit, "strategy"), "sketch")
  expect_true(!is.null(fit$betas_fixed))
})

test_that("sketch engine alias dispatches for fmri_model method", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )

  fit <- fmri_lm(
    model,
    dataset = dset,
    engine = "sketch"
  )

  expect_s3_class(fit, "fmri_lm")
  expect_equal(attr(fit, "strategy"), "sketch")
  expect_true(!is.null(fit$betas_fixed))
})

test_that("fmri_model method preserves engine-specific args for plugins", {
  engine_name <- "plugin_engine_model_test"
  captured <- NULL

  register_engine(
    engine_name,
    fit = function(model, dataset, args, cfg) {
      captured <<- list(model = model, dataset = dataset, args = args, cfg = cfg)
      Y <- as.matrix(fmridataset::get_data_matrix(dataset))
      suppressWarnings(
        fit_glm_on_transformed_series(
          model,
          Y,
          cfg = cfg,
          dataset = dataset,
          engine = engine_name,
          strategy = "engine"
        )
      )
    }
  )
  on.exit(rm(list = engine_name, envir = fmrireg:::.fmrireg_engine_registry), add = TRUE)

  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )

  fit <- fmri_lm(
    model,
    dataset = dset,
    engine = engine_name,
    engine_args = list(k = 5),
    plugin_engine_model_test = list(flag = TRUE)
  )

  expect_s3_class(fit, "fmri_lm")
  expect_true(is.list(captured))
  expect_equal(captured$args$k, 5)
  expect_true(captured$args$flag)
  expect_identical(captured$dataset, dset)
  expect_equal(attr(fit, "engine"), engine_name)
})

test_that("fmri_model method rejects unexpected arguments", {
  dset <- .demo_matrix_dataset()
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )

  expect_error(
    fmri_lm(model, dataset = dset, stray_argument = TRUE),
    "Unexpected arguments: stray_argument"
  )
})

test_that("engine capabilities reject unsupported global options before plugin fit", {
  engine_name <- "plugin_engine_caps_test"
  preflight_called <- 0L
  fit_called <- 0L

  register_engine(
    engine_name,
    preflight = function(model, dataset, args, cfg) {
      preflight_called <<- preflight_called + 1L
    },
    fit = function(model, dataset, args, cfg) {
      fit_called <<- fit_called + 1L
      fit_glm_on_transformed_series(
        model,
        as.matrix(fmridataset::get_data_matrix(dataset)),
        cfg = cfg,
        dataset = dataset,
        engine = engine_name,
        strategy = "engine"
      )
    },
    capabilities = list(
      robust = FALSE,
      preprocessing = FALSE
    )
  )
  on.exit(rm(list = engine_name, envir = fmrireg:::.fmrireg_engine_registry), add = TRUE)

  dset <- .demo_matrix_dataset()

  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = engine_name,
      robust = TRUE
    ),
    "does not support robust fitting"
  )
  expect_equal(preflight_called, 0L)
  expect_equal(fit_called, 0L)

  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = engine_name,
      volume_weights = TRUE
    ),
    "does not support volume_weights or soft_subspace preprocessing"
  )
  expect_equal(preflight_called, 0L)
  expect_equal(fit_called, 0L)
})

test_that("engine dispatcher preserves requested config and passes executed config", {
  engine_name <- "plugin_engine_exec_cfg_test"
  captured_cfg <- NULL

  register_engine(
    engine_name,
    fit = function(model, dataset, args, cfg) {
      captured_cfg <<- cfg
      fit_glm_on_transformed_series(
        model,
        as.matrix(fmridataset::get_data_matrix(dataset)),
        cfg = cfg,
        dataset = dataset,
        engine = engine_name,
        strategy = "engine"
      )
    },
    capabilities = list(
      robust = FALSE
    )
  )
  on.exit(rm(list = engine_name, envir = fmrireg:::.fmrireg_engine_registry), add = TRUE)

  dset <- .demo_matrix_dataset()
  fit <- suppressWarnings(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = engine_name,
      robust = FALSE,
      robust_options = list(max_iter = 10L)
    )
  )

  expect_s3_class(captured_cfg, "fmri_lm_config")
  expect_identical(captured_cfg$robust$type, FALSE)
  expect_identical(captured_cfg$robust$max_iter, fmri_lm_control()$robust$max_iter)

  expect_identical(attr(fit, "requested_config")$robust$max_iter, 10L)
  expect_identical(attr(fit, "executed_config")$robust$max_iter, fmri_lm_control()$robust$max_iter)
  expect_equal(attr(fit, "config"), attr(fit, "executed_config"))
})

test_that("shared engine dispatcher enforces built-in engine capabilities", {
  dset <- .demo_matrix_dataset()

  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "sketch",
      volume_weights = TRUE
    ),
    "latent_sketch does not support volume_weights or soft_subspace preprocessing"
  )

  expect_error(
    fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      engine = "rrr_gls",
      engine_args = list(rank = 1L),
      ar_voxelwise = TRUE
    ),
    "rrr_gls supports only shared \\(non-voxelwise\\) temporal covariance"
  )
})
