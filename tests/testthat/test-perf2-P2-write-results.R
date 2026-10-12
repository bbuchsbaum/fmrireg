# Performance pass 2, P2-3: write_results(format = c("h5", "nifti")) computes
# each beta / statistic volume, the mask space and each JSON sidecar once and
# hands them to every writer, instead of recomputing them per format.

.p2_wr_fit <- function() {
  set.seed(31)
  dims <- c(3, 3, 2)
  n <- 50
  scans <- lapply(1:2, function(r) {
    neuroim2::NeuroVec(array(rnorm(prod(dims) * n), c(dims, n)),
                       neuroim2::NeuroSpace(dim = c(dims, n)))
  })
  m <- array(TRUE, dims)
  m[c(2, 5, 9)] <- FALSE
  mask <- neuroim2::LogicalNeuroVol(m, neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = rep(c(5, 15, 25, 35, 45), 2),
    condition = factor(rep(c("A", "B", "C", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  dset <- neurovec_frame(scans = scans, mask = mask, TR = 1.5, event_table = etab)
  con <- contrast_set(
    pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"),
    pair_contrast(~ condition == "C", ~ condition == "B", name = "C_vs_B"),
    oneway_contrast(~ condition, name = "main")
  )
  fit <- fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~ run,
                 dataset = dset, progress = FALSE)
  list(fit = fit, mask = m)
}

.p2_f32 <- function(x) {
  readBin(writeBin(as.numeric(x), raw(), size = 4), "double", size = 4,
          n = length(x))
}

# Expected 4D volume: column j of `vals` (in-mask voxels) at volume j, 0 outside.
.p2_expected_4d <- function(vals, mask) {
  vals <- as.matrix(vals)
  out <- array(0, c(dim(mask), ncol(vals)))
  idx <- which(mask)
  for (j in seq_len(ncol(vals))) {
    sl <- array(0, dim(mask))
    sl[idx] <- vals[, j]
    out[, , , j] <- sl
  }
  out
}

.p2_h5_4d <- function(file) {
  lv <- fmristore::read_labeled_vec(file)
  labs <- names(lv)
  arr <- simplify2array(lapply(labs, function(l) as.array(lv[[l]])))
  list(labels = labs, data = array(as.vector(arr), dim(arr)))
}

.p2_nifti_4d <- function(file) {
  v <- neuroim2::read_vec(file)
  array(as.vector(as.array(v)), dim(v))
}

.p2_stat_field <- function(entry, field) {
  x <- entry[[field]]
  if (is.list(x)) x <- x[[1]]
  as.vector(x)
}

.p2_json <- function(file) {
  j <- jsonlite::read_json(file)
  j$CreationTime <- NULL
  j
}

test_that("P2-3: h5 + nifti outputs hold the fitted values (by_stat)", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  fx <- .p2_wr_fit()
  td <- tempfile("p2wr-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- write_results(fx$fit, path = td, subject = "01", task = "t",
                       space = "MNI", format = c("h5", "nifti"),
                       strategy = "by_stat",
                       contrast_stats = c("beta", "tstat", "pval", "se"))

  # Betas: every regressor, baseline included.
  B <- coef(fx$fit, type = "betas", include_baseline = TRUE)
  exp_b <- .p2_expected_4d(B, fx$mask)
  h5b <- .p2_h5_4d(res$betas$h5)
  expect_identical(h5b$labels, colnames(B))
  expect_identical(h5b$data, exp_b)
  expect_identical(.p2_nifti_4d(res$betas$nifti), .p2_f32(exp_b) |> array(dim(exp_b)))
  jb <- .p2_json(res$betas$json)
  expect_identical(unlist(jb$RegressorOrder), colnames(B))
  expect_identical(unlist(jb$DataInfo$FileFormat), c("HDF5", "NIfTI"))

  ct <- fx$fit$result$contrasts
  tcons <- ct[ct$type == "contrast", , drop = FALSE]
  fcons <- ct[ct$type == "Fcontrast", , drop = FALSE]
  expect_gt(nrow(tcons), 1L)
  expect_gt(nrow(fcons), 0L)
  checks <- list(
    list(key = "beta", rows = tcons, field = "estimate"),
    list(key = "tstat", rows = tcons, field = "stat"),
    list(key = "pval", rows = tcons, field = "prob"),
    list(key = "se", rows = tcons, field = "se"),
    list(key = "fstat", rows = fcons, field = "stat"),
    list(key = "fpval", rows = fcons, field = "prob")
  )
  for (ck in checks) {
    vals <- vapply(seq_len(nrow(ck$rows)), function(i) {
      .p2_stat_field(ck$rows$data[[i]], ck$field)
    }, numeric(sum(fx$mask)))
    exp_s <- .p2_expected_4d(vals, fx$mask)
    files <- res[[ck$key]]
    expect_false(is.null(files), label = ck$key)
    h5s <- .p2_h5_4d(files$h5)
    expect_identical(h5s$labels, ck$rows$name, label = ck$key)
    expect_identical(h5s$data, exp_s, label = ck$key)
    expect_identical(.p2_nifti_4d(files$nifti), array(.p2_f32(exp_s), dim(exp_s)),
                     label = ck$key)
    js <- .p2_json(files$json)
    expect_identical(unlist(js$ContrastOrder), ck$rows$name, label = ck$key)
    expect_identical(unlist(js$DataInfo$FileFormat), c("HDF5", "NIfTI"))
  }
})

test_that("P2-3: h5 + nifti outputs hold the fitted values (by_contrast)", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  fx <- .p2_wr_fit()
  td <- tempfile("p2wr-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  res <- write_results(fx$fit, path = td, subject = "01", task = "t",
                       space = "MNI", format = c("h5", "nifti"),
                       strategy = "by_contrast", save_betas = FALSE,
                       contrast_stats = c("beta", "tstat", "se"))
  ct <- fx$fit$result$contrasts
  for (i in seq_len(nrow(ct))) {
    nm <- ct$name[i]
    is_f <- ct$type[i] == "Fcontrast"
    stats <- if (is_f) "fstat" else c("beta", "tstat", "se")
    fields <- if (is_f) "stat" else c("estimate", "stat", "se")
    vals <- vapply(fields, function(f) .p2_stat_field(ct$data[[i]], f),
                   numeric(sum(fx$mask)))
    exp_s <- .p2_expected_4d(vals, fx$mask)
    files <- res[[nm]]
    expect_false(is.null(files), label = nm)
    h5s <- .p2_h5_4d(files$h5)
    expect_identical(h5s$labels, stats, label = nm)
    expect_identical(h5s$data, exp_s, label = nm)
    expect_identical(.p2_nifti_4d(files$nifti), array(.p2_f32(exp_s), dim(exp_s)),
                     label = nm)
    js <- .p2_json(files$json)
    expect_identical(unlist(js$StatisticOrder), stats, label = nm)
  }
})

test_that("P2-3 guard: volumes, mask space and sidecars are computed once", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("jsonlite")
  fx <- .p2_wr_fit()
  ns <- asNamespace("fmrireg")
  counts <- list()
  counting <- function(name) {
    real <- get(name, envir = ns)
    force(real)
    function(...) {
      counts[[name]] <<- (if (is.null(counts[[name]])) 0L else counts[[name]]) + 1L
      real(...)
    }
  }
  testthat::local_mocked_bindings(
    .compute_beta_volumes = counting(".compute_beta_volumes"),
    .compute_statistical_volumes = counting(".compute_statistical_volumes"),
    .fmri_dataset_mask_space = counting(".fmri_dataset_mask_space"),
    .save_betas_json_metadata = counting(".save_betas_json_metadata"),
    .save_contrasts_json_metadata = counting(".save_contrasts_json_metadata"),
    .save_single_contrast_json_metadata = counting(".save_single_contrast_json_metadata")
  )
  get_n <- function(name) if (is.null(counts[[name]])) 0L else counts[[name]]

  td <- tempfile("p2wr-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  write_results(fx$fit, path = td, subject = "01", task = "t", space = "MNI",
                format = c("h5", "nifti"), strategy = "by_stat",
                contrast_stats = c("beta", "tstat", "pval", "se"))
  # 4 t-family stats + 2 F-family stats (fstat, fpval).
  expect_identical(get_n(".compute_beta_volumes"), 1L)
  expect_identical(get_n(".compute_statistical_volumes"), 6L)
  expect_identical(get_n(".save_betas_json_metadata"), 1L)
  expect_identical(get_n(".save_contrasts_json_metadata"), 6L)
  # One mask space for the betas, one for the contrasts.
  expect_identical(get_n(".fmri_dataset_mask_space"), 2L)

  counts <- list()
  td2 <- tempfile("p2wr-")
  dir.create(td2)
  on.exit(unlink(td2, recursive = TRUE), add = TRUE)
  write_results(fx$fit, path = td2, subject = "01", task = "t", space = "MNI",
                format = c("h5", "nifti"), strategy = "by_contrast",
                save_betas = FALSE, contrast_stats = c("beta", "tstat"))
  expect_identical(get_n(".save_single_contrast_json_metadata"),
                   nrow(fx$fit$result$contrasts))
  expect_identical(get_n(".fmri_dataset_mask_space"), 1L)
})
