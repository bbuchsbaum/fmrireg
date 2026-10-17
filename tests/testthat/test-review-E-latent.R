# Review 2026-10, stream E2: fmri_latent_lm must honour (or refuse) the
# requested noise model instead of silently fitting OLS.

.review_e_latent_data <- function(seed = 3L) {
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  rl <- 120L; n <- 2L * rl; r <- 4L; V <- 60L
  ons <- seq(10, 220, 18)
  ev <- data.frame(
    onset = c(ons, ons),
    condition = factor(rep(c("a", "b"), length.out = 2L * length(ons))),
    run = rep(1:2, each = length(ons))
  )
  Z <- matrix(0, n, r)
  for (j in seq_len(r)) Z[, j] <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  L <- matrix(rnorm(V * r), V, r)
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(V, 1, 1)),
                                    neuroim2::NeuroSpace(c(V, 1, 1)))
  list(
    Z = Z, L = L, ev = ev, rl = rl,
    latent = latent_frame(list(basis = Z, loadings = L, mask = mask), TR = 2,
                          run_length = c(rl, rl), event_table = ev),
    dense = matrix_frame(Z %*% t(L), TR = 2, run_length = c(rl, rl), event_table = ev)
  )
}

test_that("E2: latent AR(1) fit differs from OLS and matches the dense AR fit", {
  d <- .review_e_latent_data()
  fit_ols <- fmri_latent_lm(onset ~ hrf(condition), block = ~run,
                            dataset = d$latent, durations = 0)
  fit_ar <- fmri_latent_lm(onset ~ hrf(condition), block = ~run,
                           dataset = d$latent, durations = 0, autocor = "ar1")

  B_ols <- as.matrix(coef(fit_ols))
  B_ar <- as.matrix(coef(fit_ar))
  expect_gt(max(abs(B_ar - B_ols)), 1e-3)

  # Dense reference: the same joint GLS solve on the reconstructed voxel data
  # with run-pooled AR(1) (one phi per run shared across voxels), estimated
  # from the joint-design OLS residuals as the latent fit does. (The built-in
  # fmri_lm() backends estimate phi from per-run OLS residuals instead.)
  X <- as.matrix(design_matrix(fit_ar$model))
  runs <- split(seq_len(nrow(X)), rep(1:2, each = d$rl))
  fit_dense <- fmrireg:::solve_integrated_glm(
    X, d$Z %*% t(d$L),
    fmri_lm_control(estimation = estimation_spec("joint"), noise = noise_spec("ar1")),
    run_indices = runs
  )
  phi_dense <- unlist(fit_dense$ar_coef)
  phi_lat <- unlist(fit_ar$result$ar_coef)
  expect_length(phi_lat, 2L)
  expect_equal(unname(phi_lat), unname(phi_dense), tolerance = 1e-6)

  # Voxel betas reconstructed from the component betas equal the dense betas.
  B_comp <- t(fit_ar$result$betas$data[[1]]$estimate[[1]])          # p x r
  expect_equal(unname(B_comp %*% t(d$L)), unname(fit_dense$betas[, , drop = FALSE]),
               tolerance = 1e-6)

  # Component-level SEs match an explicit GLS fit on whitened component scores.
  w <- fmrireg:::ar_whiten_transform(X, d$Z, phi = as.list(phi_lat), run_indices = runs)
  ref <- stats::lm.fit(w$X, w$Y)
  sig2 <- colSums(ref$residuals^2) / ref$df.residual
  se_ref <- sqrt(outer(diag(chol2inv(qr.R(ref$qr)))[order(ref$qr$pivot)], sig2))
  se_lat <- t(fit_ar$result$betas$data[[1]]$se[[1]])
  expect_equal(unname(se_lat), unname(se_ref), tolerance = 1e-6)
})

test_that("E2: robust latent fits error instead of silently fitting OLS", {
  d <- .review_e_latent_data()
  expect_error(
    fmri_latent_lm(onset ~ hrf(condition), block = ~run, dataset = d$latent,
                   durations = 0, robust = TRUE),
    "not supported for latent"
  )
  expect_error(
    fmri_lm(onset ~ hrf(condition), block = ~run, dataset = d$latent,
            control = fmri_lm_control(estimation = estimation_spec("joint"),
                                      robust = robust_spec("huber")),
            compute = compute_spec(voxel_chunks = 1L)),
    "not supported for latent"
  )
})

test_that("E2: unsupported latent configurations error", {
  d <- .review_e_latent_data()
  # bootstrap was accepted and silently ignored
  expect_error(
    fmri_latent_lm(onset ~ hrf(condition), block = ~run, dataset = d$latent,
                   durations = 0, bootstrap = TRUE),
    "not supported for latent"
  )
  # Internal guard (MA and voxelwise AR are already rejected upstream for
  # joint fits; the latent fitter checks independently).
  cfg <- fmri_lm_control(estimation = estimation_spec("joint"),
                         noise = noise_spec("ar1"))
  cfg_ma <- cfg; cfg_ma$ar$q <- 1L
  expect_error(fmrireg:::.latent_check_supported_config(cfg_ma), "not supported for latent")
  cfg_vx <- cfg; cfg_vx$ar$voxelwise <- TRUE
  expect_error(fmrireg:::.latent_check_supported_config(cfg_vx), "not supported for latent")
  expect_true(fmrireg:::.latent_check_supported_config(cfg))
})

test_that("E2: latent AR(2) also matches the dense fit", {
  d <- .review_e_latent_data(seed = 9L)
  fit_ar <- fmri_latent_lm(onset ~ hrf(condition), block = ~run,
                           dataset = d$latent, durations = 0, autocor = "ar2")
  X <- as.matrix(design_matrix(fit_ar$model))
  fit_dense <- fmrireg:::solve_integrated_glm(
    X, d$Z %*% t(d$L),
    fmri_lm_control(estimation = estimation_spec("joint"), noise = noise_spec("ar2")),
    run_indices = split(seq_len(nrow(X)), rep(1:2, each = d$rl))
  )
  expect_equal(unlist(fit_ar$result$ar_coef), unlist(fit_dense$ar_coef),
               tolerance = 1e-6, ignore_attr = TRUE)
  B_comp <- t(fit_ar$result$betas$data[[1]]$estimate[[1]])
  expect_equal(unname(B_comp %*% t(d$L)), unname(fit_dense$betas[, , drop = FALSE]),
               tolerance = 1e-6)
})
