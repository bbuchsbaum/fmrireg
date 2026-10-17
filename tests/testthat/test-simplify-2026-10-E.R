# Regression tests for simplification stream E (2026-10): print methods,
# with_package(), BIDS export, imports.

.e_two_run_frame <- function() {
  set.seed(11)
  n <- 60
  ev <- data.frame(
    onset = rep(c(5, 20, 35, 50), 2),
    condition = factor(rep(c("A", "B", "C", "A"), 2)),
    run = rep(1:2, each = 4)
  )
  matrix_frame(matrix(rnorm(2 * n * 6), 2 * n, 6), TR = 1,
               run_length = c(n, n), event_table = ev)
}

test_that("print.fmri_lm reports the fit's scope, stored df, formula and contrasts (B2)", {
  dset <- .e_two_run_frame()
  cons <- list(
    fmridesign::pair_contrast(~ condition == "A", ~ condition == "B", name = "AvB"),
    fmridesign::oneway_contrast(~ condition, name = "main_cond")
  )
  fit <- fmri_lm(onset ~ hrf(condition, contrasts = cons), block = ~run,
                 dataset = dset)
  out <- capture.output(print(fit))

  df_stored <- unique(fit[["result"]][["df"]][["inference"]])
  expect_length(df_stored, 1L)
  expect_true(any(grepl("Estimation scope: runwise_meta", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Inference df (residual): %s", format(df_stored)),
                        out, fixed = TRUE)))
  expect_true(any(grepl("onset ~ hrf(condition, contrasts = cons)", out, fixed = TRUE)))
  expect_true(any(grepl("Simple contrasts (1): AvB", out, fixed = TRUE)))
  expect_true(any(grepl("F-contrasts (1): main_cond", out, fixed = TRUE)))
  expect_true(any(grepl("Parameters: 3 event + 8 baseline", out, fixed = TRUE)))

  # joint (chunkwise) scope, AR and robust configuration lines
  fit_ar <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                    strategy = "chunkwise", nchunks = 2, cor_struct = "ar1")
  out_ar <- capture.output(print(fit_ar))
  expect_true(any(grepl("Estimation scope: joint", out_ar, fixed = TRUE)))
  expect_true(any(grepl("AR structure: ar1", out_ar, fixed = TRUE)))
  expect_true(any(grepl(
    sprintf("Inference df (residual): %s", format(unique(fit_ar[["result"]][["df"]][["inference"]]))),
    out_ar, fixed = TRUE)))

  fit_rob <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                     robust = "huber")
  out_rob <- capture.output(print(fit_rob))
  expect_true(any(grepl("Robust method: huber (tuning constant 1.345)", out_rob, fixed = TRUE)))

  # fmri_rlm keeps its banner and falls through to the fmri_lm summary
  fit_rlm <- fmri_rlm(onset ~ hrf(condition), block = ~run, dataset = dset)
  out_rlm <- capture.output(print(fit_rlm))
  expect_true(any(grepl("Robust fMRI Linear Model", out_rlm, fixed = TRUE)))
  expect_true(any(grepl("Estimation scope: runwise_meta", out_rlm, fixed = TRUE)))
})

test_that("print.fmri_lm names the engine and its df for engine fits (B2)", {
  dset <- .e_two_run_frame()
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 engine = "latent_sketch")
  out <- capture.output(print(fit))
  df_stored <- unique(fit[["result"]][["df"]][["inference"]])
  expect_true(any(grepl("Engine: latent_sketch", out, fixed = TRUE)))
  expect_true(any(grepl(
    sprintf("Inference df (%s): %s", fit[["result"]][["df"]][["method"]],
            format(df_stored, digits = 5L)),
    out, fixed = TRUE)))
})

test_that("print.fmri_model shows event terms, event count, baseline basis and total columns (B3)", {
  dset <- .e_two_run_frame()
  base <- baseline_model(basis = "poly", degree = 2,
                         sframe = fmrihrf::sampling_frame(c(60, 60), TR = 1))
  fm <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                          baseline_model = base, dataset = dset)
  out <- capture.output(print(fm))

  expect_true(any(grepl("Formula: onset ~ hrf(condition)", out, fixed = TRUE)))
  expect_true(any(grepl("condition (<event_term>)", out, fixed = TRUE)))
  expect_true(any(grepl("Number of Events: 8", out, fixed = TRUE)))
  expect_true(any(grepl("Basis type: poly", out, fixed = TRUE)))
  expect_true(any(grepl("Degree: 2", out, fixed = TRUE)))
  # 3 condition columns + (2 poly + 1 intercept) per run x 2 runs
  expect_identical(ncol(design_matrix(fm)), 9L)
  expect_true(any(grepl("Total design columns: 9", out, fixed = TRUE)))
})

test_that("with_package checks every package and names all missing ones (B6)", {
  expect_null(fmrireg:::with_package(c("stats", "utils")))
  err <- tryCatch(
    fmrireg:::with_package(c("stats", "notapkg.zzz1", "notapkg.zzz2")),
    error = function(e) conditionMessage(e)
  )
  expect_identical(
    err,
    "Please install the notapkg.zzz1, notapkg.zzz2 packages to use this functionality"
  )
  # single missing package keeps the original wording
  expect_error(fmrireg:::with_package("notapkg.zzz1"),
               "Please install the notapkg.zzz1 package to use this functionality",
               fixed = TRUE)
})

# ---- write_results(): pinned file layout and sidecars ------------------------

.e_bids_fit <- function() {
  set.seed(123)
  dims <- c(3, 3, 2); nt <- 50
  scans <- lapply(1:2, function(run) {
    arr <- array(rnorm(prod(dims) * nt), c(dims, nt))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, nt)))
  })
  mask_arr <- array(TRUE, dims)
  mask_arr[1, 1, 1] <- FALSE
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim = dims))
  ev <- data.frame(
    onset = rep(c(5, 15, 25, 35, 45), 2),
    condition = factor(rep(c("A", "B", "C", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  dset <- neurovec_frame(scans = scans, mask = mask, TR = 1.5, event_table = ev)
  con <- contrast_set(
    pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"),
    oneway_contrast(~ condition, name = "cond_F")
  )
  fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~run, dataset = dset,
          durations = 0, strategy = "chunkwise", nchunks = 2, progress = FALSE)
}

.e_expected_files <- function(stems, fmt) {
  exts <- c(h5 = ".h5", nifti = ".nii.gz")[fmt]
  sort(c(paste0("sub-01_task-t_space-MNI_", stems, ".json"),
         as.vector(outer(paste0("sub-01_task-t_space-MNI_", stems), exts, paste0))))
}

test_that("write_results writes the pinned file set for every format/layout/desc", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("RNifti")
  fit <- .e_bids_fit()
  all_stats <- c("beta", "tstat", "pval", "se")

  stems <- list(
    "by_stat:GLM:all" = c("desc-GLM_betas", "desc-beta_bold", "desc-fpval_bold",
                          "desc-fstat_bold", "desc-pval_bold", "desc-se_bold",
                          "desc-tstat_bold"),
    "by_stat:run1:all" = c("desc-run1_betas", "desc-run1beta_bold",
                           "desc-run1fpval_bold", "desc-run1fstat_bold",
                           "desc-run1pval_bold", "desc-run1se_bold",
                           "desc-run1tstat_bold"),
    "by_contrast:GLM:all" = c("contrast-AvsB_desc-GLM_bold",
                              "contrast-condF_desc-GLM_bold", "desc-GLM_betas"),
    "by_contrast:run1:all" = c("contrast-AvsB_desc-run1_bold",
                               "contrast-condF_desc-run1_bold", "desc-run1_betas"),
    "by_stat:GLM:tstat" = c("desc-GLM_betas", "desc-fstat_bold", "desc-tstat_bold"),
    "by_stat:run1:tstat" = c("desc-run1_betas", "desc-run1fstat_bold",
                             "desc-run1tstat_bold"),
    "by_contrast:GLM:tstat" = c("contrast-AvsB_desc-tstat_bold",
                                "contrast-condF_desc-fstat_bold", "desc-GLM_betas"),
    "by_contrast:run1:tstat" = c("contrast-AvsB_desc-run1tstat_bold",
                                 "contrast-condF_desc-run1fstat_bold",
                                 "desc-run1_betas")
  )

  for (key in names(stems)) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1]]
    for (fmt in list("h5", "nifti", c("h5", "nifti"))) {
      out <- tempfile("e_bids_")
      dir.create(out)
      res <- suppressWarnings(write_results(
        fit, path = out, subject = "01", task = "t", space = "MNI",
        desc = parts[2], format = fmt, strategy = parts[1],
        contrast_stats = if (parts[3] == "all") all_stats else "tstat"
      ))
      expect_identical(sort(list.files(out, recursive = TRUE, all.files = TRUE)),
                       .e_expected_files(stems[[key]], fmt),
                       info = paste(key, paste(fmt, collapse = "+")))
      # every returned path exists and nothing written is unreported
      returned <- unlist(res, use.names = FALSE)
      expect_setequal(normalizePath(returned),
                      normalizePath(list.files(out, full.names = TRUE)))
      unlink(out, recursive = TRUE)
    }
  }
})

test_that("write_results sidecars and maps carry the pinned contents", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("RNifti")
  fit <- .e_bids_fit()
  out <- tempfile("e_bids_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  res <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                       format = c("h5", "nifti"))
  read <- function(f) jsonlite::read_json(file.path(out, paste0("sub-01_task-t_space-MNI_", f)))

  # returned structure: betas, then one entry per statistic, files in
  # h5/json/nifti order
  expect_identical(names(res), c("betas", "beta", "tstat", "pval", "se", "fstat", "fpval"))
  expect_identical(names(res$betas), c("h5", "json", "nifti"))
  expect_identical(names(res$tstat), c("h5", "json", "nifti"))

  bj <- read("desc-GLM_betas.json")
  beta_mat <- coef(fit, type = "betas", include_baseline = TRUE)
  expect_identical(unlist(bj$RegressorOrder), colnames(beta_mat))
  expect_identical(bj$ModelInfo$NumRegressors, ncol(beta_mat))
  expect_identical(unlist(bj$DataInfo$FileFormat), c("HDF5", "NIfTI"))
  expect_identical(unlist(bj$DataInfo$Storage), c("LabeledVolumeSet", "NIfTI-1 image"))
  expect_identical(bj$Description, "Raw regressor beta coefficients from GLM analysis")
  expect_identical(bj$SpatialReference, "MNI")

  tj <- read("desc-tstat_bold.json")
  expect_identical(tj$StatisticType, "TSTAT")
  expect_identical(unlist(tj$ContrastOrder), "A_vs_B")
  expect_identical(tj$ContrastDefinitions$A_vs_B$Type, "t-contrast")
  expect_identical(tj$DataInfo$Units, "t-statistic")
  expect_identical(tj$Description, "Contrast tstat statistics from GLM analysis")
  expect_equal(tj$DegreesOfFreedom, unique(fit$result$df$inference))

  fj <- read("desc-fstat_bold.json")
  expect_identical(fj$StatisticType, "FSTAT")
  expect_identical(fj$ContrastDefinitions$cond_F$Type, "F-contrast")
  expect_identical(fj$DataInfo$Units, "F-statistic")

  # NIfTI map values equal the fitted statistics (float32) inside the mask
  mask <- as.vector(fmrireg:::.fmri_dataset_mask_space(fit$dataset)$mask_array)
  tvol <- neuroim2::read_vec(res$tstat$nifti)
  tstat <- fit$result$contrasts$data[[which(fit$result$contrasts$name == "A_vs_B")]]$stat
  expect_equal(as.array(tvol)[, , , 1][mask], tstat, tolerance = 1e-6)
  expect_true(all(as.array(tvol)[, , , 1][!mask] == 0))
  bvol <- neuroim2::read_vec(res$betas$nifti)
  expect_equal(matrix(as.array(bvol), ncol = ncol(beta_mat))[mask, ],
               unname(beta_mat), tolerance = 1e-6, ignore_attr = TRUE)

  # by_contrast sidecar for a multi-statistic contrast
  out2 <- tempfile("e_bids_")
  dir.create(out2)
  on.exit(unlink(out2, recursive = TRUE), add = TRUE)
  write_results(fit, path = out2, subject = "01", task = "t", space = "MNI",
                format = "nifti", strategy = "by_contrast", save_betas = FALSE)
  cj <- jsonlite::read_json(file.path(out2, "sub-01_task-t_space-MNI_contrast-AvsB_desc-GLM_bold.json"))
  expect_identical(unlist(cj$StatisticOrder), c("beta", "tstat", "pval", "se"))
  expect_identical(cj$ContrastInfo$Type, "t-contrast")
  expect_identical(cj$DataInfo$FileFormat, "NIfTI")
  fcj <- jsonlite::read_json(file.path(out2, "sub-01_task-t_space-MNI_contrast-condF_desc-GLM_bold.json"))
  expect_identical(unlist(fcj$StatisticOrder), c("fstat", "fpval"))
  expect_identical(fcj$ContrastInfo$Type, "F-contrast")
})

test_that("write_results with h5 + nifti computes each map and writes each sidecar once (B9)", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("RNifti")
  fit <- .e_bids_fit()
  json_paths <- character(0)
  n_beta <- 0L
  n_stat <- 0L
  real_json <- jsonlite::write_json
  real_beta <- fmrireg:::.compute_beta_volumes
  real_stat <- fmrireg:::.compute_statistical_volumes
  local_mocked_bindings(
    write_json = function(x, path, ...) {
      json_paths <<- c(json_paths, basename(path))
      real_json(x, path, ...)
    },
    .package = "jsonlite"
  )
  local_mocked_bindings(
    .compute_beta_volumes = function(...) { n_beta <<- n_beta + 1L; real_beta(...) },
    .compute_statistical_volumes = function(...) { n_stat <<- n_stat + 1L; real_stat(...) }
  )

  out <- tempfile("e_bids_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                format = c("h5", "nifti"), strategy = "by_stat")

  # 1 beta sidecar + 4 t-family + 2 F-family statistic sidecars
  expect_length(json_paths, 7L)
  expect_false(anyDuplicated(json_paths) > 0)
  expect_identical(n_beta, 1L)
  expect_identical(n_stat, 6L)
})

test_that("voxel plots read the observed series of the requested voxel", {
  dset <- .e_two_run_frame()
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset)
  Y <- as.matrix(fmrireg:::.dset_data_matrix(dset))
  for (v in c(1L, 4L)) {
    parts <- fmrireg:::.voxel_fit_parts(fit, v)
    expect_identical(parts$y, as.numeric(Y[, v]))
  }
})
