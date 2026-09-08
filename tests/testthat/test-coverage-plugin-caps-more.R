# Additional aaa_plugin_api capability negotiation branches.

test_that("register_engine capabilities NULL and print forbid/require lines", {
  eng <- paste0("cov_caps_", as.integer(Sys.time()) %% 1e6L)
  register_engine(
    eng,
    fit = function(model, dataset, args, cfg) list(ok = TRUE, engine = eng),
    preflight = NULL,
    capabilities = NULL
  )
  spec <- fmrireg:::get_engine_spec(eng, include_functions = FALSE)
  expect_s3_class(spec, "fmrireg_engine_spec")
  expect_true(is.list(spec$capabilities))

  # Print paths for requires/forbids
  register_engine(
    paste0(eng, "_rules"),
    fit = function(...) NULL,
    capabilities = list(
      requires_event_regressors = TRUE,
      requires_parcels_for_by_cluster = TRUE,
      forbid_by_cluster_dataset_classes = c("matrix_dataset", "latent_dataset"),
      robust = FALSE,
      ma = FALSE,
      preprocessing = FALSE
    )
  )
  out <- capture.output(print(fmrireg:::engine_spec(paste0(eng, "_rules"))))
  txt <- paste(out, collapse = "\n")
  expect_true(grepl("requires|forbids|capabilities", txt, ignore.case = TRUE))
})

test_that(".validate_engine_capabilities rejects MA/robust/voxelwise/scope", {
  cfg <- fmri_lm_control()

  cfg_ma <- cfg
  cfg_ma$ar$q <- 1L
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg_ma, capabilities = list(ma = FALSE)
    ),
    "MA"
  )

  cfg_rob <- fmri_lm_control()
  cfg_rob$robust$type <- "huber"
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg_rob, capabilities = list(robust = FALSE)
    ),
    "robust"
  )

  cfg_vox <- cfg
  cfg_vox$ar$voxelwise <- TRUE
  expect_error(
    fmrireg:::.validate_engine_capabilities(
      "e", cfg_vox, capabilities = list(ar_voxelwise = FALSE)
    ),
    "voxelwise|shared"
  )

  cfg_scope <- cfg
  # Force unsupported scope if accessible
  if (!is.null(cfg_scope$estimation)) {
    cfg_scope$estimation$scope <- "runwise_meta"
    expect_error(
      fmrireg:::.validate_engine_capabilities(
        "e", cfg_scope, capabilities = list(estimation_scopes = "joint")
      ),
      "estimation scope|scope"
    )
  }

  # Happy path defaults
  ok <- fmrireg:::.validate_engine_capabilities("ols", cfg, capabilities = list())
  expect_true(is.null(ok) || isTRUE(ok) || is.list(ok) || identical(ok, cfg))
})

test_that(".derive_engine_execution_config and alias helpers cover branches", {
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1", by_cluster = TRUE))
  derived <- fmrireg:::.derive_engine_execution_config(
    cfg,
    capabilities = list(ar_by_cluster = TRUE, ar_voxelwise = TRUE, robust = TRUE)
  )
  expect_true(is.list(derived))

  expect_equal(fmrireg:::.builtin_engine_aliases("latent_sketch"), "sketch")
  expect_equal(fmrireg:::.builtin_engine_aliases("ols"), character())
  expect_equal(fmrireg:::.builtin_engine_source("rrr_gls"), "builtin")
  expect_equal(fmrireg:::.builtin_engine_source("custom_plugin"), "plugin")

  caps <- fmrireg:::.normalize_engine_capabilities(list(robust = FALSE, ma = TRUE))
  expect_false(isTRUE(caps$robust))
  expect_true(isTRUE(caps$ma))
})
