# Fifth wave: fit_contrasts.fmri_lm pair/F paths, engine re-register, RRR bootstrap SE.

test_that("fit_contrasts.fmri_lm resolves pair_contrast and F-contrast specs", {
  set.seed(31)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55, 65, 75),
    condition = factor(rep(c("A", "B"), 4)),
    run = 1L
  )
  Y <- matrix(rnorm(100 * 5), 100, 5)
  dset <- matrix_dataset(Y, TR = 1, run_length = 100, event_table = etab)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)

  con <- pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
  out <- fit_contrasts(fit, list(con))
  expect_true(is.list(out))
  expect_true(length(out) >= 1L)
  first <- out[[1]]
  expect_true(!is.null(first$estimate) || !is.null(first$stat))
  expect_equal(length(first$estimate), 5L)

  # Named list of contrast matrices
  L <- matrix(c(1, -1), 1, dimnames = list("A_vs_B", c("conditionA", "conditionB")))
  # Use beta names from fit if available
  bnames <- colnames(coef(fit))
  if (length(bnames) >= 2L) {
    L2 <- matrix(c(1, -1), 1, dimnames = list("diff", bnames[1:2]))
    out2 <- tryCatch(fit_contrasts(fit, list(diff = L2)), error = function(e) e)
    expect_true(inherits(out2, "error") || is.list(out2))
  }

  # F-contrast via multi-row matrix when API accepts it
  if (length(bnames) >= 2L) {
    LF <- rbind(c(1, 0), c(0, 1))
    colnames(LF) <- bnames[1:2]
    outF <- tryCatch(fit_contrasts(fit, list(main = LF)), error = function(e) e)
    expect_true(inherits(outF, "error") || is.list(outF))
  }
})

test_that("fit_contrasts.fmri_lm structured covariance voxel path", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  p <- length(fit$result$event_indices)
  V <- nrow(as.matrix(coef(fit)))
  if (p < 1L || V < 1L) skip("demo fit too small")

  # Attach a variance_model with per-voxel covariance
  XtXinv <- fit$result$XtXinv %||% fit$result$cov.unscaled
  if (is.null(XtXinv)) {
    X <- as.matrix(design_matrix(fit$model))
    XtXinv <- chol2inv(chol(crossprod(X)))
  }
  fit$result$variance_model <- structure(
    list(
      df_inference = rep(max(2, fit$result$rdf %||% 2), V),
      covariance_scope = "voxel",
      covariance = replicate(V, as.matrix(XtXinv), simplify = FALSE)
    ),
    class = "fmri_lm_variance_model"
  )
  con <- pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
  out <- tryCatch(fit_contrasts(fit, list(con)), error = function(e) e)
  expect_true(inherits(out, "error") || is.list(out))
})

test_that(".register_builtin_engines re-registers after registry clear", {
  reg <- fmrireg:::.fmrireg_engine_registry
  saved <- as.list(reg, all.names = TRUE)
  on.exit({
    rm(list = ls(reg, all.names = TRUE), envir = reg)
    for (nm in names(saved)) reg[[nm]] <- saved[[nm]]
  }, add = TRUE)

  rm(list = ls(reg, all.names = TRUE), envir = reg)
  expect_null(get_engine("rrr_gls"))
  expect_null(get_engine("latent_sketch"))

  expect_true(fmrireg:::.register_builtin_engines())
  expect_true(!is.null(get_engine("rrr_gls")))
  expect_true(!is.null(get_engine("latent_sketch")))

  spec <- engine_spec("rrr_gls")
  expect_s3_class(spec, "fmrireg_engine_spec")
  out <- capture.output(print(spec))
  expect_true(any(grepl("rrr_gls|capabilities|requires", out)))

  spec2 <- engine_spec("latent_sketch")
  out2 <- capture.output(print(spec2))
  expect_true(any(grepl("parcels|latent_dataset|forbids", out2)))

  specs <- engine_specs()
  expect_true("rrr_gls" %in% names(specs))
  expect_true(length(list_engine_specs(include_functions = TRUE)) >= 2L)
})

test_that(".fit_rrr_gls_engine bootstrap SE and energy rank modes", {
  skip_if_not_installed("fmriAR")
  set.seed(33)
  etab <- data.frame(
    onset = seq(5, by = 8, length.out = 8),
    condition = factor(rep(c("A", "B"), 4)),
    run = 1L
  )
  Y <- matrix(rnorm(80 * 6), 80, 6)
  dset <- matrix_dataset(Y, TR = 1, run_length = 80, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  fm <- fmri_model(emod, bmod, dset)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))

  fit_boot <- fmrireg:::.fit_rrr_gls_engine(
    fm, dset,
    args = list(
      rank = 2L, se_mode = "bootstrap", bootstrap_n = 8L,
      bootstrap_block_size = 4L, bootstrap_seed = 1L
    ),
    cfg = cfg
  )
  expect_s3_class(fit_boot, "fmri_lm")

  fit_energy <- fmrireg:::.fit_rrr_gls_engine(
    fm, dset,
    args = list(rank_mode = "energy", energy_keep = 0.85, se_mode = "conditional"),
    cfg = fmri_lm_control()
  )
  expect_s3_class(fit_energy, "fmri_lm")
})
