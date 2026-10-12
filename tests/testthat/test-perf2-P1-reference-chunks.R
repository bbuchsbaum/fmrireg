# Performance pass 2, item P1-3: the reference ("slow") chunkwise backend
# estimates the noise model (AR plan, robust weights) once over all voxels and
# hands it to every voxel chunk, so results do not depend on `voxel_chunks`,
# and the AR estimation / IRLS run once rather than once per chunk.

.p13_fixture <- function(V = 90L, n_run = 70L, seed = 13L) {
  set.seed(seed)
  onsets <- seq(4, 2 * n_run - 20, by = 16)
  ev <- data.frame(
    onset = c(onsets, onsets),
    cond = factor(rep(c("A", "B"), length.out = 2 * length(onsets))),
    run = rep(1:2, each = length(onsets))
  )
  Y <- vapply(seq_len(V), function(v) {
    c(as.numeric(arima.sim(list(ar = 0.5), n_run)),
      as.numeric(arima.sim(list(ar = 0.3), n_run))) + 10
  }, numeric(2 * n_run))
  Y[c(9, 50, 101), ] <- Y[c(9, 50, 101), ] + 6   # outlier volumes
  list(dset = matrix_frame(Y, TR = 2, run_length = c(n_run, n_run),
                           event_table = ev), Y = Y)
}

.p13_fit <- function(fx, chunks, noise = noise_spec("ar1"),
                     robust = robust_spec("none")) {
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AvB"))
  fmri_lm(onset ~ hrf(cond, contrasts = con), block = ~run, dataset = fx$dset,
          durations = 0,
          control = fmri_lm_control(estimation = estimation_spec("joint"),
                                    noise = noise, robust = robust),
          compute = compute_spec(voxel_chunks = chunks, backend = "reference"))
}

.p13_summary <- function(fit) {
  b <- fit$result$betas$data[[1]]
  con <- fit$result$contrasts$data[[1]]
  list(estimate = b$estimate[[1]], se = b$se[[1]], stat = b$stat[[1]],
       con_est = con$estimate, con_se = con$se, con_stat = con$stat,
       sigma = as.numeric(fit$result[["sigma"]]))
}

.p13_cases <- list(
  ar1 = list(noise = noise_spec("ar1"), robust = robust_spec("none")),
  robust_run = list(noise = noise_spec("iid"),
                    robust = robust_spec("huber", scale_scope = "run")),
  robust_voxel = list(noise = noise_spec("iid"),
                      robust = robust_spec("bisquare", scale_scope = "voxel",
                                           max_iter = 4L)),
  ar1_robust = list(noise = noise_spec("ar1"),
                    robust = robust_spec("huber", scale_scope = "run"))
)

test_that("reference backend results do not depend on voxel_chunks", {
  fx <- .p13_fixture()
  for (nm in names(.p13_cases)) {
    cs <- .p13_cases[[nm]]
    ref <- .p13_summary(.p13_fit(fx, 1L, cs$noise, cs$robust))
    for (k in c(3L, 7L)) {
      got <- .p13_summary(.p13_fit(fx, k, cs$noise, cs$robust))
      for (fld in names(ref)) {
        expect_equal(got[[fld]], ref[[fld]], tolerance = 1e-10,
                     info = paste(nm, k, fld))
      }
    }
  }
})

test_that("reference backend AR fit equals GLS with the single estimated plan", {
  # With the noise model fixed, each voxel's fit is GLS: lm on the whitened
  # data, with the whitening from the one AR plan of the joint fit.
  fx <- .p13_fixture(V = 12L)
  fit <- .p13_fit(fx, 1L)
  X <- as.matrix(design_matrix(fit$model))
  phi <- fit$result$ar_coef
  phi <- if (is.list(phi[[1]])) phi[[1]] else phi
  runs <- list(1:70, 71:140)
  wh <- fmrireg:::ar_whiten_transform(X, fx$Y, phi, FALSE, run_indices = runs)
  b <- fit$result$betas$data[[1]]$estimate[[1]]
  se <- fit$result$betas$data[[1]]$se[[1]]
  for (v in c(1L, 7L, 12L)) {
    ref <- lm(wh$Y[, v] ~ wh$X - 1)
    keep <- seq_len(ncol(b))
    expect_equal(unname(b[v, ]), unname(coef(ref))[keep], tolerance = 1e-8)
    expect_equal(unname(se[v, ]),
                 unname(summary(ref)$coefficients[keep, 2]), tolerance = 1e-8)
  }
})

test_that("reference backend estimates AR and robust weights once per fit", {
  fx <- .p13_fixture(V = 30L)
  ar_calls <- 0L
  irls_calls <- 0L
  orig_est <- fmrireg:::.estimate_ar_via_fmriAR
  orig_fitter <- fmrireg:::robust_iterative_fitter
  testthat::local_mocked_bindings(
    .estimate_ar_via_fmriAR = function(...) {
      ar_calls <<- ar_calls + 1L
      orig_est(...)
    },
    robust_iterative_fitter = function(...) {
      if (is.null(list(...)$weights_fixed)) irls_calls <<- irls_calls + 1L
      orig_fitter(...)
    }
  )
  .p13_fit(fx, 5L)
  expect_equal(ar_calls, 1L)
  ar_calls <- 0L
  .p13_fit(fx, 5L, noise = noise_spec("ar1"),
           robust = robust_spec("huber", scale_scope = "run"))
  expect_equal(ar_calls, 1L)
  expect_equal(irls_calls, 1L)
})
