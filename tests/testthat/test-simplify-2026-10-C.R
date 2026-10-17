# Regression tests for the 2026-10 simplification, stream C (engines, plugin
# API, latent fits, betas).

# --- B7: rrr_gls resolves censoring like every other fitting path ---------

# Two runs whose censored volumes hold gross artefacts, so the AR estimate
# moves visibly when censoring is (or is not) applied.
.rrr_censor_fixture <- function(bad = c(10L, 11L, 12L, 55L, 56L)) {
  set.seed(42)
  runs <- c(40L, 40L)
  total <- sum(runs)
  datamat <- matrix(rnorm(total * 6L), total, 6L)
  datamat[bad, ] <- datamat[bad, ] + 50
  ev <- do.call(rbind, lapply(seq_along(runs), function(r) {
    on <- seq(4, runs[r] * 2 - 12, by = 8)
    data.frame(onset = on,
               condition = factor(rep(c("A", "B"), length.out = length(on))),
               run = r)
  }))
  cen <- rep(FALSE, total)
  cen[bad] <- TRUE
  list(
    censored = matrix_frame(datamat, TR = 2, run_length = runs,
                            event_table = ev, censor = cen),
    plain = matrix_frame(datamat, TR = 2, run_length = runs, event_table = ev),
    bad = bad
  )
}

.rrr_censor_fit <- function(ds, censor = NULL) {
  fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = ds,
    engine = "rrr_gls", engine_args = list(rank = 1L),
    control = fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec("ar1", censor = censor)
    )
  )
}

.rrr_phi <- function(fit) unlist(fit[["ar_coef"]], use.names = FALSE)
.rrr_coef <- function(fit) as.numeric(fit[["betas_fixed"]])

test_that("rrr_gls uses a frame-carried censor column (B7)", {
  fx <- .rrr_censor_fixture()

  from_column <- .rrr_censor_fit(fx$censored)
  explicit <- .rrr_censor_fit(fx$plain, censor = fx$bad)
  auto <- .rrr_censor_fit(fx$censored, censor = "auto")
  uncensored <- .rrr_censor_fit(fx$plain)

  # The censor column is honoured exactly as the explicit indices are.
  expect_equal(.rrr_phi(from_column), .rrr_phi(explicit), tolerance = 1e-12)
  expect_equal(.rrr_coef(from_column), .rrr_coef(explicit), tolerance = 1e-12)
  expect_equal(.rrr_phi(auto), .rrr_phi(explicit), tolerance = 1e-12)
  expect_equal(.rrr_coef(auto), .rrr_coef(explicit), tolerance = 1e-12)
  # ... and censoring the artefacts changes the AR estimate.
  expect_false(isTRUE(all.equal(.rrr_phi(explicit), .rrr_phi(uncensored))))
})

test_that("rrr_gls censor = 'none' ignores the frame column (B7)", {
  fx <- .rrr_censor_fixture()

  # Before the fix the raw string reached fmriAR as as.integer("none").
  expect_no_warning(
    opted_out <- .rrr_censor_fit(fx$censored, censor = "none")
  )
  uncensored <- .rrr_censor_fit(fx$plain)

  expect_equal(.rrr_phi(opted_out), .rrr_phi(uncensored), tolerance = 1e-12)
  expect_equal(.rrr_coef(opted_out), .rrr_coef(uncensored), tolerance = 1e-12)
})

# --- B8: standard_error(<fmri_latent_lm>, recon = TRUE) after an AR fit ------

.latent_se_fixture <- function(autocor) {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("neuroim2")
  set.seed(42)
  n_time <- 60L
  n_comp <- 4L
  n_vox <- 20L
  basis <- matrix(rnorm(n_time * n_comp), n_time, n_comp)
  for (t in 2:n_time) basis[t, ] <- basis[t, ] + 0.5 * basis[t - 1L, ]
  loadings <- matrix(rnorm(n_vox * n_comp), n_vox, n_comp)
  lvec <- suppressWarnings(fmristore::LatentNeuroVec(
    basis = basis, loadings = loadings,
    space = neuroim2::NeuroSpace(c(5, 4, 1, n_time)),
    mask = rep(TRUE, n_vox), offset = rep(0, n_vox)
  ))
  lds <- latent_frame(lvec, TR = 1, run_length = n_time, event_table = data.frame(
    onset = c(8, 20, 32, 44), condition = factor(c("A", "B", "A", "B")), run = 1L
  ))
  fit <- fmri_latent_lm(
    onset ~ hrf(condition, contrasts = pair_contrast(
      ~ condition == "A", ~ condition == "B", name = "AvB")),
    block = ~ run, dataset = lds, durations = 0, autocor = autocor
  )
  list(fit = fit, Z = basis, L = loadings,
       X = as.matrix(design_matrix(fit$model)))
}

# Closed form: with whitened design Xw and scores Zw, Var(B[i, ]) =
# C[i, i] * Sigma, Sigma = Rw'Rw / df the component residual covariance, so a
# voxel with loadings l has se = sqrt(C[i, i] * l' Sigma l).
.latent_se_closed_form <- function(Xw, Zw, L, cmat = NULL) {
  C <- solve(crossprod(Xw))
  B <- C %*% crossprod(Xw, Zw)
  Rw <- Zw - Xw %*% B
  Sigma <- crossprod(Rw) / (nrow(Xw) - ncol(Xw))
  q <- rowSums((L %*% Sigma) * L)
  if (is.null(cmat)) {
    sapply(1:2, function(i) sqrt(C[i, i] * q))
  } else {
    sqrt(drop(t(cmat) %*% C %*% cmat) * q)
  }
}

test_that("latent AR fit reconstructs standard errors in voxel space (B8)", {
  fx <- .latent_se_fixture("ar1")
  phi <- unlist(fx$fit$result[["ar_coef"]])
  expect_length(phi, 1L)

  # AR(1) whitening without exact_first: row 1 kept, rows t >= 2 get
  # x_t - phi * x_{t-1}.
  W <- diag(nrow(fx$X))
  W[cbind(2:nrow(W), 1:(nrow(W) - 1L))] <- -phi
  Xw <- W %*% fx$X
  Zw <- W %*% fx$Z

  se <- as.matrix(standard_error(fx$fit, type = "betas", recon = TRUE))
  expect_equal(unname(se), .latent_se_closed_form(Xw, Zw, fx$L),
               tolerance = 1e-8)

  se_con <- as.matrix(standard_error(fx$fit, type = "contrasts", recon = TRUE))
  expect_equal(unname(se_con[, 1]),
               .latent_se_closed_form(Xw, Zw, fx$L, cmat = c(1, -1, 0, 0, 0, 0)[1:ncol(fx$X)]),
               tolerance = 1e-8)
})

test_that("latent IID standard errors in voxel space are unchanged (B8)", {
  fx <- .latent_se_fixture("none")

  se <- as.matrix(standard_error(fx$fit, type = "betas", recon = TRUE))
  expect_equal(unname(se), .latent_se_closed_form(fx$X, fx$Z, fx$L),
               tolerance = 1e-8)

  # The formula used before the fix: cov() of the OLS residuals rescaled to
  # the residual df (the design has an intercept, so the centring is a no-op).
  R <- fx$Z - fx$X %*% qr.coef(qr(fx$X), fx$Z)
  CR <- cov(R) * (nrow(R) - 1) / (nrow(R) - ncol(fx$X))
  C <- solve(crossprod(fx$X))
  old <- sapply(1:2, function(i) sqrt(rowSums((fx$L %*% (CR * C[i, i])) * fx$L)))
  expect_equal(unname(se), old, tolerance = 1e-10)
})

# --- glm_lss reports non-finite betas accurately --------------------------

test_that("glm_lss names non-finite trial betas instead of a Cholesky failure", {
  ev <- data.frame(onset = c(10, 30), condition = factor(c("A", "B")), run = 1)
  dset <- matrix_frame(matrix(rnorm(100), 50, 2), TR = 2, run_length = 50,
                       event_table = ev)
  emod <- event_model(onset ~ hrf(condition), data = ev, block = ~ run,
                      sampling_frame = fmrihrf::sampling_frame(50, TR = 2))
  local_mocked_bindings(
    estimate_betas = function(...) list(betas_ran = matrix(c(1, NaN, 2, 3), 2, 2))
  )
  expect_error(glm_lss(dset, emod, fmrihrf::HRF_SPMG1),
               "non-finite trial betas \\(1 of 4 values\\)")
})
