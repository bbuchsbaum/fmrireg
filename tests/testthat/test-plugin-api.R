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
