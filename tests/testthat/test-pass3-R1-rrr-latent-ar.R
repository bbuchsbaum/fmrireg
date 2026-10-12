# Pass 3, R1-2: rrr_gls on latent data with run/global-pooled AR noise stays
# in latent space. Voxel residual cross-products are Z_res (L'L) Z_res', so
# the pooled AR estimate can use Z_res Q Lambda^(1/2) (L'L = Q Lambda Q');
# whitening is linear in the response, so the whitened fit stays latent too.
# The result must equal the fit of the dense reconstruction Z L'.

.p3r_latent_fixture <- function(n_time = 100L, nrun = 2L, n_comp = 5L,
                                n_vox = 60L, phi = c(0.5, -0.1), seed = 41L) {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  ons <- seq(6, 2 * n_time - 20, by = 14)
  ev <- do.call(rbind, lapply(seq_len(nrun), function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("a", "b", "c"), length.out = length(ons))),
               run = r)
  }))
  ntot <- n_time * nrun
  basis <- matrix(0, ntot, n_comp)
  for (r in seq_len(nrun)) for (j in seq_len(n_comp)) {
    basis[(r - 1L) * n_time + seq_len(n_time), j] <-
      as.numeric(arima.sim(list(ar = phi), n_time))
  }
  basis[, 1] <- basis[, 1] + 2
  # Correlated, unequally scaled loadings: L'L is far from the identity.
  loadings <- matrix(stats::rnorm(n_vox * n_comp), n_vox, n_comp) %*%
    matrix(c(1, 0.6, 0, 0, 0,
             0, 2, 0.5, 0, 0,
             0, 0, 0.5, 0.3, 0,
             0, 0, 0, 3, 1,
             0, 0, 0, 0, 1), n_comp, n_comp)
  lvec <- suppressWarnings(fmristore::LatentNeuroVec(
    basis = basis, loadings = loadings,
    space = neuroim2::NeuroSpace(c(n_vox, 1, 1, ntot)),
    mask = rep(TRUE, n_vox), offset = rep(0, n_vox)
  ))
  lds <- latent_frame(lvec, TR = 2, run_length = rep(n_time, nrun), event_table = ev)
  dense <- matrix_frame(basis %*% t(loadings), TR = 2,
                        run_length = rep(n_time, nrun), event_table = ev)
  list(lds = lds, dense = dense, n_vox = n_vox)
}

.p3r_fit <- function(dset, noise, engine_args = list(rank = 2L)) {
  fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
          engine = "rrr_gls", engine_args = engine_args,
          control = fmri_lm_control(noise = noise))
}

.p3r_same <- function(f_lat, f_den, tol = 1e-8) {
  expect_equal(f_lat$result$ar_coef, f_den$result$ar_coef, tolerance = tol)
  expect_equal(as.matrix(coef(f_lat)), as.matrix(coef(f_den)), tolerance = tol)
  expect_equal(as.matrix(standard_error(f_lat)), as.matrix(standard_error(f_den)),
               tolerance = tol)
  expect_equal(as.matrix(stats(f_lat)), as.matrix(stats(f_den)), tolerance = tol)
  expect_equal(f_lat$betas_fixed, f_den$betas_fixed, tolerance = tol)
  expect_equal(f_lat$sigma2, f_den$sigma2, tolerance = tol)
  expect_equal(f_lat$result$rss, f_den$result$rss, tolerance = tol)
  expect_equal(f_lat$rrr$singular_values, f_den$rrr$singular_values, tolerance = tol)
}

test_that("latent rrr_gls with pooled AR equals the dense fit (phi, betas, SEs)", {
  fx <- .p3r_latent_fixture()
  noises <- list(
    noise_spec("ar1", pooling = "run"),
    noise_spec("ar1", pooling = "global"),
    noise_spec("ar2", pooling = "run"),
    noise_spec("ar1", pooling = "run", shared_estimator = "mean_series"),
    noise_spec("ar1", pooling = "global", exact_first = TRUE)
  )
  for (nz in noises) {
    f_lat <- .p3r_fit(fx$lds, nz)
    f_den <- .p3r_fit(fx$dense, nz)
    expect_gt(max(abs(unlist(f_den$result$ar_coef))), 0.2)
    .p3r_same(f_lat, f_den)
  }
  # Full rank (= GLS) as well.
  nz <- noise_spec("ar1", pooling = "run")
  .p3r_same(.p3r_fit(fx$lds, nz, list()), .p3r_fit(fx$dense, nz, list()))
})

test_that("guard: latent rrr_gls with pooled AR never reconstructs the dense T x V response", {
  fx <- .p3r_latent_fixture(n_vox = 30L)
  local_mocked_bindings(
    .rrr_extract_response_matrix = function(dataset) {
      stop("dense T x V reconstruction of the response")
    }
  )
  for (nz in list(noise_spec("ar1", pooling = "run"),
                  noise_spec("ar2", pooling = "global"),
                  noise_spec("ar1", shared_estimator = "mean_series"))) {
    expect_no_error(f <- .p3r_fit(fx$lds, nz))
    expect_equal(ncol(f$betas_fixed), fx$n_vox)
  }
})
