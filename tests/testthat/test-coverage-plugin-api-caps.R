# aaa_plugin_api.R: registration, capability context, basis inject, GLM helpers.

make_plugin_model <- function(n = 50L, V = 4L, seed = 31) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(5, 15, 25, 35),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "constant", sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    Y = Y,
    cfg = fmri_lm_control()
  )
}

test_that("register_engine validation and engine_spec listing/print", {
  expect_error(register_engine("", function(...) NULL), "nzchar|length")
  expect_error(register_engine("bad_fit", fit = "nope"), "must be a function")
  expect_error(
    register_engine("bad_pre", fit = function(...) NULL, preflight = 1),
    "NULL or a function"
  )

  eng_name <- paste0("cov_eng_", as.integer(Sys.time()) %% 100000L)
  register_engine(
    eng_name,
    fit = function(model, dataset, args, cfg) list(ok = TRUE),
    preflight = function(...) TRUE,
    capabilities = list(robust = FALSE, ma = TRUE, ar_by_cluster = FALSE)
  )
  expect_true(!is.null(get_engine(eng_name)))
  expect_null(get_engine(""))
  expect_null(get_engine(1))

  spec <- get_engine_spec(eng_name, include_functions = TRUE)
  expect_s3_class(spec, "fmrireg_engine_spec")
  expect_true(is.function(spec$fit))
  expect_true(is.function(spec$preflight))
  expect_false(spec$capabilities$robust)

  expect_null(get_engine_spec(""))
  expect_null(get_engine_spec("definitely_missing_engine_xyz"))

  specs <- list_engine_specs()
  expect_true(eng_name %in% names(specs))
  expect_true("rrr_gls" %in% names(engine_specs()))

  out <- capture.output(print(engine_spec(eng_name)))
  expect_true(any(grepl(eng_name, out)))
  expect_true(any(grepl("capabilities", out)))

  # Builtin latent_sketch print shows parcels requirement
  ls_out <- capture.output(print(engine_spec("latent_sketch")))
  expect_true(any(grepl("parcels|forbids|requires", ls_out, ignore.case = TRUE)))
})

test_that(".validate_engine_context covers parcels and forbid-class rules", {
  fx <- make_plugin_model()
  caps <- list(
    requires_event_regressors = TRUE,
    requires_parcels_for_by_cluster = TRUE,
    forbid_by_cluster_dataset_classes = "matrix_dataset"
  )
  cfg <- fx$cfg
  cfg$ar$by_cluster <- TRUE

  expect_error(
    fmrireg:::.validate_engine_context(
      "latent_sketch", fx$model, fx$dataset, args = list(lowrank = list()),
      cfg = cfg, capabilities = caps
    ),
    "parcels"
  )

  expect_error(
    fmrireg:::.validate_engine_context(
      "latent_sketch", fx$model, fx$dataset,
      args = list(lowrank = list(parcels = 1:4)),
      cfg = cfg, capabilities = caps
    ),
    "by_cluster|matrix_dataset"
  )

  # Happy path without by_cluster
  cfg2 <- fx$cfg
  ok <- fmrireg:::.validate_engine_context(
    "ols", fx$model, fx$dataset, args = list(), cfg = cfg2,
    capabilities = list(requires_event_regressors = TRUE)
  )
  expect_true(is.list(ok))
})

test_that(".validate_engine_capabilities rejects variance/df/scope/preprocessing", {
  cfg <- fmri_lm_control()
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg, capabilities = list(variance_methods = "sandwich")
    ),
    "variance method"
  )
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg, capabilities = list(df_methods = "satterthwaite")
    ),
    "df method"
  )
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg, capabilities = list(estimation_scopes = "runwise_meta")
    ),
    "estimation scope"
  )

  cfg_prep <- fmri_lm_control()
  cfg_prep$volume_weights$enabled <- TRUE
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg_prep, capabilities = list(preprocessing = FALSE)
    ),
    "volume_weights|preprocessing"
  )

  cfg_cluster <- fmri_lm_control()
  cfg_cluster$ar$by_cluster <- TRUE
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg_cluster, capabilities = list(ar_by_cluster = FALSE)
    ),
    "parcel-specific|by_cluster|AR"
  )

  derived <- fmrireg:::.derive_engine_execution_config(
    cfg_cluster,
    capabilities = list(ar_by_cluster = FALSE, ar_voxelwise = FALSE, robust = FALSE)
  )
  expect_false(isTRUE(derived$ar$by_cluster) || isTRUE(derived$noise$by_cluster))
})

test_that(".fmrireg_inject_registered_bases / replace_basis_expr rewrite formulas", {
  register_basis("cov_basis_inject", function(...) fmrihrf::HRF_SPMG1)

  f <- onset ~ hrf(condition, basis = "cov_basis_inject")
  injected <- fmrireg:::.fmrireg_inject_registered_bases(f)
  expect_true(inherits(injected, "formula"))
  # Injected RHS should mention resolve_basis
  rhs <- deparse(injected[[3]])
  expect_true(any(grepl("resolve_basis", rhs)))

  # Non-formula passthrough
  expect_equal(fmrireg:::.fmrireg_inject_registered_bases(1:3), 1:3)

  # Unregistered basis left alone
  f2 <- onset ~ hrf(condition, basis = "not_registered_xyz")
  left <- fmrireg:::.fmrireg_inject_registered_bases(f2)
  expect_true(grepl("not_registered_xyz", paste(deparse(left), collapse = "")))
})

test_that("fit_glm_on_transformed_series / from_suffstats / external dims", {
  fx <- make_plugin_model(n = 40L, V = 3L)
  fit <- fit_glm_on_transformed_series(fx$model, fx$Y, cfg = fx$cfg)
  expect_s3_class(fit, "fmri_lm")
  expect_equal(ncol(coef(fit, type = "betas")), 3L)

  expect_error(
    fit_glm_on_transformed_series(list(), fx$Y),
    "fmri_model"
  )
  expect_error(
    fit_glm_on_transformed_series(fx$model, fx$Y, cfg = list()),
    "fmri_lm_control"
  )
  cfg_ar <- fmri_lm_control(ar_options = list(struct = "ar1"))
  expect_error(
    fit_glm_on_transformed_series(fx$model, fx$Y, cfg = cfg_ar),
    "IID OLS only"
  )

  X <- as.matrix(design_matrix(fx$model))
  XtX <- crossprod(X)
  XtS <- crossprod(X, fx$Y)
  StS <- colSums(fx$Y^2)
  suff <- fit_glm_from_suffstats(
    fx$model, XtX, XtS, StS, df = nrow(X) - ncol(X)
  )
  expect_s3_class(suff, "fmri_lm")

  expect_null(fmrireg:::.external_response_dataset_dims(NULL))
  dims <- fmrireg:::.external_response_dataset_dims(fx$dataset)
  expect_equal(dims, c(40L, 3L))

  expect_true(fmrireg:::.validate_external_response_dataset(NULL, c(40, 3), "ctx"))
  expect_error(
    fmrireg:::.validate_external_response_dataset(fx$dataset, c(10, 3), "ctx"),
    "dimensions"
  )
})

test_that("fit_glm_with_config covers AR/robust on tiny Y", {
  skip_if_not_installed("fmriAR")
  fx <- make_plugin_model(n = 60L, V = 3L, seed = 41)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  fit <- fit_glm_with_config(fx$model, fx$Y, cfg = cfg)
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))
})
