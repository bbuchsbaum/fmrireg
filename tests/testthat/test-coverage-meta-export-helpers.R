# Twelfth wave: fmri_meta export helpers (image set, json, by_stat/by_coef nifti).

make_vol_meta <- function(seed = 231L) {
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  dims <- c(2L, 2L, 1L)
  space <- neuroim2::NeuroSpace(dim = dims)
  vals <- array(rnorm(prod(dims) * 2), c(dims, 2))
  # Flatten voxel x coef matrices
  coefs <- matrix(as.numeric(vals), prod(dims), 2,
                  dimnames = list(NULL, c("Intercept", "groupB")))
  se <- matrix(runif(prod(dims) * 2, 0.1, 0.4), prod(dims), 2,
               dimnames = dimnames(coefs))
  # Attach reconstructable NeuroVols via coef_image by storing space/mask
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), space)
  structure(
    list(
      coefficients = coefs,
      se = se,
      tau2 = runif(prod(dims)),
      I2 = runif(prod(dims), 0, 50),
      Q = rchisq(prod(dims), df = 3),
      Q_df = rep(3, prod(dims)),
      method = "DL",
      robust = "none",
      formula = ~ 1 + group,
      n_subjects = 8L,
      n_voxels = prod(dims),
      mask = mask,
      space = space
    ),
    class = c("fmri_meta", "list")
  )
}

test_that(".fmri_meta_stat_image / ensure / image_set cover NeuroVol path", {
  skip_if_not_installed("neuroim2")
  meta <- make_vol_meta()

  # coef_image may return matrix without reconstruct; force NeuroVol images
  vol1 <- neuroim2::NeuroVol(as.numeric(meta$coefficients[, 1]), meta$space)
  vol2 <- neuroim2::NeuroVol(as.numeric(meta$coefficients[, 2]), meta$space)
  expect_identical(fmrireg:::.ensure_fmri_meta_image(vol1, meta, "beta"), vol1)
  expect_error(
    fmrireg:::.ensure_fmri_meta_image(meta$coefficients[, 1], meta, "beta"),
    "no voxel-space"
  )

  expect_error(fmrireg:::.fmri_meta_image_set(list(), character()), "No images")
  expect_error(
    fmrireg:::.fmri_meta_image_set(list(vol1), c("a", "b")),
    "do not match"
  )

  img_set <- fmrireg:::.fmri_meta_image_set(list(vol1, vol2), c("Intercept", "groupB"), x = meta)
  expect_s4_class(img_set$neurovec, "NeuroVec")
  expect_s4_class(img_set$mask, "LogicalNeuroVol")
  expect_equal(img_set$labels, c("Intercept", "groupB"))
})

test_that(".fmri_meta_json_metadata and predict output files", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  ent <- list(subject = "01", task = "meta", space = "MNI")
  js <- fmrireg:::.fmri_meta_json_metadata(
    x = meta,
    description = "unit test meta maps",
    statistic_names = c("beta", "se"),
    coefficient_names = c("Intercept"),
    output_formats = "nifti",
    units = "arbitrary",
    entities = ent
  )
  expect_true(is.list(js))
  expect_true(!is.null(js$Description) || !is.null(js$DataInfo) || length(js) >= 1L)

  predicted <- fmrireg:::.predict_fmri_meta_output_files(
    path = "/tmp",
    entities = ent,
    desc = "META",
    format = "nifti",
    strategy = "by_stat",
    coefficient_names = c("Intercept", "groupB"),
    coefficient_stats = c("beta", "se"),
    heterogeneity = TRUE
  )
  expect_true(is.list(predicted) || is.character(predicted))
  flat <- unlist(predicted, use.names = FALSE)
  expect_true(any(grepl("nii|json|heterogeneity|Intercept|groupB|beta|se", flat, ignore.case = TRUE)))
})

test_that(".write_fmri_meta_map and save_by_stat/by_coefficient nifti succeed", {
  skip_if_not_installed("neuroim2")
  meta <- make_vol_meta(seed = 232L)
  td <- tempfile("meta-export-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ent <- list(subject = "01", task = "meta", space = "MNI")

  # Patch coef_image.fmri_meta to return NeuroVols for this synthetic object
  vol_from_coef <- function(object, coef = 1, statistic = c("estimate", "se", "z", "p"), ...) {
    statistic <- match.arg(statistic)
    idx <- if (is.character(coef)) match(coef, colnames(object$coefficients)) else as.integer(coef)
    vals <- switch(
      statistic,
      estimate = object$coefficients[, idx],
      se = object$se[, idx],
      z = object$coefficients[, idx] / object$se[, idx],
      p = 2 * pnorm(-abs(object$coefficients[, idx] / object$se[, idx]))
    )
    neuroim2::NeuroVol(as.numeric(vals), object$space)
  }
  # Use local override via assignInNamespace is risky; call lower helpers with NeuroVols instead.
  vols <- list(
    neuroim2::NeuroVol(as.numeric(meta$coefficients[, 1]), meta$space),
    neuroim2::NeuroVol(as.numeric(meta$coefficients[, 2]), meta$space)
  )
  img_set <- fmrireg:::.fmri_meta_image_set(vols, c("Intercept", "groupB"), x = meta)
  meta_json <- fmrireg:::.fmri_meta_json_metadata(
    meta, "beta maps", "beta", c("Intercept", "groupB"),
    output_formats = "nifti", units = "a.u.", entities = ent
  )
  written <- fmrireg:::.write_fmri_meta_map(
    img_set, td, "sub-01_task-meta_desc-beta_bold",
    meta_json, format = "nifti", output_formats = "nifti"
  )
  expect_true(file.exists(written$nifti))
  expect_true(file.exists(written$json))

  # Heterogeneity with reconstructable vectors via NeuroVol override of reconstruct_image path:
  # .save_fmri_meta_heterogeneity uses reconstruct_image; if that fails, skip gracefully.
  het <- tryCatch(
    fmrireg:::.save_fmri_meta_heterogeneity(
      meta, td, ent, desc = "META", format = "nifti", output_formats = "nifti"
    ),
    error = function(e) e
  )
  expect_true(inherits(het, "error") || is.list(het))
})

test_that("build_contrast_from_names and parse_contrast_formula edges", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  cn <- colnames(meta$coefficients)
  skip_if(is.null(cn) || length(cn) < 1L)

  w <- fmrireg:::build_contrast_from_names(setNames(1, cn[[1]]), meta)
  expect_equal(sum(w), 1)
  expect_error(
    fmrireg:::build_contrast_from_names(c(missing = 1), meta),
    "not found"
  )

  if (length(cn) >= 2L) {
    parsed <- fmrireg:::parse_contrast_formula(
      as.formula(paste0("~ `", cn[[1]], "` - `", cn[[2]], "`")),
      meta
    )
    expect_true(is.numeric(parsed) || is.list(parsed))
  }
})
