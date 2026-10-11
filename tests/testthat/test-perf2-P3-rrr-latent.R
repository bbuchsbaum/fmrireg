# Performance pass 2, stream P3-4: rrr_gls on a latent dataset stays in
# latent space (scores Z, loadings L) instead of forming the dense T x V
# response Z L'. Every step of the iid fit is linear in the response, so the
# result must equal the fit of the dense reconstruction.

p3_latent_fixture <- function(n_time = 120L, n_comp = 6L, n_vox = 40L, seed = 21L) {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  ev <- data.frame(
    onset = seq(6, 220, by = 14),
    condition = factor(rep(c("a", "b", "c"), length.out = length(seq(6, 220, by = 14)))),
    run = 1L
  )
  basis <- matrix(stats::rnorm(n_time * n_comp), n_time, n_comp)
  basis[, 1] <- basis[, 1] + 2
  loadings <- matrix(stats::rnorm(n_vox * n_comp), n_vox, n_comp)
  lvec <- fmristore::LatentNeuroVec(
    basis = basis, loadings = loadings,
    space = neuroim2::NeuroSpace(c(n_vox, 1, 1, n_time)),
    mask = rep(TRUE, n_vox), offset = rep(0, n_vox)
  )
  lds <- latent_frame(lvec, TR = 2, run_length = n_time, event_table = ev)
  Y <- basis %*% t(loadings)
  dense <- matrix_frame(Y, TR = 2, run_length = n_time, event_table = ev)
  list(lds = lds, dense = dense, Y = Y, ev = ev)
}

p3_rrr_fit <- function(dset, engine_args = list(), ...) {
  fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset,
          engine = "rrr_gls", engine_args = engine_args, ...)
}

p3_expect_same_fit <- function(f_lat, f_den, tol = 1e-8) {
  expect_equal(as.matrix(coef(f_lat)), as.matrix(coef(f_den)), tolerance = tol)
  expect_equal(as.matrix(standard_error(f_lat)), as.matrix(standard_error(f_den)),
               tolerance = tol)
  expect_equal(as.matrix(stats(f_lat)), as.matrix(stats(f_den)), tolerance = tol)
  expect_equal(f_lat$betas_fixed, f_den$betas_fixed, tolerance = tol)
  expect_equal(f_lat$sigma2, f_den$sigma2, tolerance = tol)
  expect_equal(f_lat$result$rss, f_den$result$rss, tolerance = tol)
  expect_equal(f_lat$rrr$singular_values, f_den$rrr$singular_values, tolerance = tol)
  expect_equal(f_lat$rrr$rank_used, f_den$rrr$rank_used)
}

test_that("P3-4: latent rrr_gls equals the dense reconstruction (full and reduced rank)", {
  fx <- p3_latent_fixture()
  cases <- list(
    list(),
    list(rank = 1L),
    list(rank = 2L),
    list(rank_mode = "energy", energy_keep = 0.6)
  )
  for (a in cases) {
    f_lat <- p3_rrr_fit(fx$lds, a)
    f_den <- p3_rrr_fit(fx$dense, a)
    p3_expect_same_fit(f_lat, f_den)
  }
  # contrasts on task parameters
  con <- fmri_lm(onset ~ hrf(condition, contrasts = pair_contrast(~ condition == "a",
                                                                   ~ condition == "b",
                                                                   name = "a_b")),
                 block = ~ run, dataset = fx$lds, engine = "rrr_gls",
                 engine_args = list(rank = 2L))
  con_d <- fmri_lm(onset ~ hrf(condition, contrasts = pair_contrast(~ condition == "a",
                                                                     ~ condition == "b",
                                                                     name = "a_b")),
                   block = ~ run, dataset = fx$dense, engine = "rrr_gls",
                   engine_args = list(rank = 2L))
  expect_equal(as.matrix(stats(con, type = "contrasts")),
               as.matrix(stats(con_d, type = "contrasts")), tolerance = 1e-8)
  expect_equal(as.matrix(standard_error(con, type = "contrasts")),
               as.matrix(standard_error(con_d, type = "contrasts")), tolerance = 1e-8)
})

test_that("P3-4: latent full-rank rrr_gls matches lm() on a reconstructed voxel", {
  fx <- p3_latent_fixture()
  f_lat <- p3_rrr_fit(fx$lds)
  X <- as.matrix(design_matrix(f_lat$model))
  for (v in c(1L, 17L, 40L)) {
    s <- summary(stats::lm(fx$Y[, v] ~ X - 1))$coefficients
    expect_equal(unname(as.matrix(coef(f_lat))[v, ]), unname(s[1:3, 1]), tolerance = 1e-8)
    expect_equal(unname(as.matrix(standard_error(f_lat))[v, ]), unname(s[1:3, 2]),
                 tolerance = 1e-8)
  }
})

test_that("P3-4: latent bootstrap SEs equal the dense bootstrap with the same seed", {
  fx <- p3_latent_fixture()
  a <- list(rank = 2L, se_mode = "bootstrap", bootstrap_n = 25L, bootstrap_seed = 9L)
  p3_expect_same_fit(p3_rrr_fit(fx$lds, a), p3_rrr_fit(fx$dense, a))
})

test_that("P3-4: latent AR rrr_gls still equals the dense fit", {
  fx <- p3_latent_fixture()
  f_lat <- p3_rrr_fit(fx$lds, list(rank = 2L), ar_options = list(struct = "ar1"))
  f_den <- p3_rrr_fit(fx$dense, list(rank = 2L), ar_options = list(struct = "ar1"))
  p3_expect_same_fit(f_lat, f_den)
})

test_that("P3-4 guard: latent iid rrr_gls never reconstructs the dense T x V response", {
  fx <- p3_latent_fixture()
  local_mocked_bindings(
    .rrr_extract_response_matrix = function(dataset) {
      stop("dense T x V reconstruction of the response")
    }
  )
  expect_no_error(f1 <- p3_rrr_fit(fx$lds, list(rank = 2L)))
  expect_no_error(f2 <- p3_rrr_fit(
    fx$lds, list(rank = 2L, se_mode = "bootstrap", bootstrap_n = 5L, bootstrap_seed = 1L)
  ))
  expect_equal(dim(f1$betas_fixed), c(ncol(design_matrix(f1$model)), ncol(fx$Y)))
})
