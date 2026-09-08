# Latent lm reconstruction paths (coef/SE/stats with recon=TRUE).

.make_latent_fit <- function() {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(42)
  n_time <- 60
  n_comp <- 4
  n_voxels <- 20
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  loadings <- matrix(rnorm(n_voxels * n_comp), n_voxels, n_comp)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = neuroim2::NeuroSpace(c(5, 4, 1, n_time)),
    mask = rep(TRUE, n_voxels),
    offset = rep(0, n_voxels)
  )
  lds <- fmridataset::latent_dataset(
    source = list(lvec), TR = 1, run_length = n_time
  )
  lds$event_table <- data.frame(
    onset = c(8, 20, 32, 44),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )

  fit <- fmri_latent_lm(
    onset ~ hrf(condition), block = ~ run, dataset = lds,
    durations = 0, autocor = "none", robust = FALSE
  )
  list(fit = fit, n_voxels = n_voxels, n_comp = n_comp)
}

test_that("fmri_latent_lm reconstructs coef/se/stats into voxel space", {
  fx <- .make_latent_fit()
  fit <- fx$fit

  cf <- coef(fit, type = "estimates", recon = TRUE)
  expect_s3_class(cf, "tbl_df")
  expect_equal(nrow(cf), fx$n_voxels)
  expect_true(all(is.finite(as.matrix(cf))))

  cf_betas <- coef(fit, type = "betas", recon = TRUE)
  expect_equal(as.matrix(cf_betas), as.matrix(cf))

  # Two conditions and four components: subsetting to two components used
  # to match n_coefficients == n_selected and take the legacy row branch,
  # yielding voxels × 4 instead of voxels × 2.
  b_latent <- as.matrix(coef(fit, type = "estimates", recon = FALSE))
  expect_equal(nrow(b_latent), 2L)
  expect_equal(ncol(b_latent), fx$n_comp)

  lvec <- if (!is.null(fit$dataset$lvec)) {
    fit$dataset$lvec
  } else {
    fit$dataset$backend$data[[1]]
  }
  lds <- lvec@loadings
  expect_equal(as.matrix(cf), lds %*% t(b_latent), tolerance = 1e-8)
  expect_equal(ncol(cf), 2L)

  cf_sub <- coef(fit, type = "estimates", recon = TRUE, comp = 1:2)
  expect_equal(nrow(cf_sub), fx$n_voxels)
  expect_equal(ncol(cf_sub), 2L)
  expect_equal(as.matrix(cf_sub), lds[, 1:2, drop = FALSE] %*% t(b_latent[, 1:2, drop = FALSE]),
               tolerance = 1e-8)

  expect_error(
    coef(fit, type = "estimates", recon = TRUE, comp = fx$n_comp + 1L),
    regexp = "."
  )

  se <- standard_error(fit, type = "estimates", recon = TRUE)
  expect_s3_class(se, "tbl_df")
  expect_equal(nrow(se), fx$n_voxels)
  expect_true(all(as.matrix(se) > 0, na.rm = TRUE))

  st <- stats(fit, type = "estimates", recon = TRUE)
  expect_s3_class(st, "tbl_df")
  expect_equal(dim(st), dim(cf))
  expect_equal(as.matrix(st), as.matrix(cf) / as.matrix(se), tolerance = 1e-8)

  se_c <- standard_error(fit, type = "contrasts", recon = TRUE)
  expect_equal(nrow(se_c), 0L)
})

test_that("fmri_latent_lm reconstructs contrast SE via stubbed contrasts", {
  fx <- .make_latent_fit()
  fit <- fx$fit
  fit2 <- fit
  fit2$result$contrasts <- tibble::tibble(
    type = "contrast",
    name = "A_vs_B",
    conmat = list(matrix(c(1, -1), ncol = 1)),
    colind = list(fit$result$event_indices[1:2]),
    data = list(list())
  )
  se_forced <- standard_error(fit2, type = "contrasts", recon = TRUE)
  expect_equal(nrow(se_forced), fx$n_voxels)
  expect_true("A_vs_B" %in% names(se_forced))
  expect_true(all(is.finite(as.matrix(se_forced))))
})

test_that("latent recon keeps type-based orientation when n_comp equals n_conditions", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")

  set.seed(7)
  n_time <- 80
  n_comp <- 4
  n_voxels <- 16
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  loadings <- matrix(rnorm(n_voxels * n_comp), n_voxels, n_comp)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = neuroim2::NeuroSpace(c(4, 4, 1, n_time)),
    mask = rep(TRUE, n_voxels),
    offset = rep(0, n_voxels)
  )
  lds <- fmridataset::latent_dataset(
    source = list(lvec), TR = 1, run_length = n_time
  )
  lds$event_table <- data.frame(
    onset = c(8, 20, 32, 44, 56, 68),
    condition = factor(c("A", "B", "C", "D", "A", "B")),
    run = 1L
  )
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"))
  fit <- fmri_latent_lm(
    onset ~ hrf(condition, contrasts = con), block = ~ run, dataset = lds,
    durations = 0, autocor = "none", robust = FALSE
  )

  b_latent <- as.matrix(coef(fit, type = "estimates", recon = FALSE))
  expect_equal(dim(b_latent), c(4L, n_comp))
  cf <- coef(fit, type = "estimates", recon = TRUE)
  expect_equal(dim(as.matrix(cf)), c(n_voxels, 4L))
  expect_equal(as.matrix(cf), loadings %*% t(b_latent), tolerance = 1e-8)

  c_latent <- as.matrix(coef(fit, type = "contrasts", recon = FALSE))
  expect_equal(nrow(c_latent), n_comp)
  cf_c <- coef(fit, type = "contrasts", recon = TRUE)
  expect_equal(nrow(cf_c), n_voxels)
  expect_equal(as.matrix(cf_c), loadings %*% c_latent, tolerance = 1e-8)
})

test_that("coef.fmri_latent_lm errors when LatentNeuroVec is missing", {
  fx <- .make_latent_fit()
  fit <- fx$fit
  fit$dataset$lvec <- NULL
  if (!is.null(fit$dataset$backend)) {
    fit$dataset$backend$data <- NULL
  }
  expect_error(
    coef(fit, type = "estimates", recon = TRUE),
    "Cannot find LatentNeuroVec"
  )
})
