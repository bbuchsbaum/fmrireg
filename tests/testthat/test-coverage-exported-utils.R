# Meaningful coverage for previously untested exported utilities,
# boundary cases, and error paths.

tiny_matrix_dataset <- function(n_time = 60L, n_vox = 8L, n_events = 6L, TR = 1) {
  set.seed(42)
  onsets <- seq(5, n_time - 10, length.out = n_events)
  etab <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = n_events)),
    run = 1L
  )
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  list(
    dataset = matrix_dataset(Y, TR = TR, run_length = n_time, event_table = etab),
    events = etab
  )
}

test_that("generate_main_effect_contrast and interaction contrasts have expected dims", {
  des <- expand.grid(
    Time = factor(1:4),
    Cond = factor(c("face", "scene"))
  )

  main <- generate_main_effect_contrast(des, "Time")
  expect_equal(dim(main), c(8L, 3L))
  expect_true(all(abs(colSums(main)) < 1e-10))

  inter <- generate_interaction_contrast(des, c("Time", "Cond"))
  expect_equal(dim(inter), c(8L, 3L))

  expect_error(
    generate_main_effect_contrast(des, c("Time", "Cond")),
    "exactly one factor"
  )
  expect_error(
    generate_interaction_contrast(des, "Missing"),
    "factors %in% names"
  )
})

test_that("hrf_smoothing_kernel returns normalized gram/cosine kernels", {
  form <- onset ~ trialwise(basis = "gaussian")
  K_gram <- hrf_smoothing_kernel(40, TR = 1.5, form = form, method = "gram")
  expect_equal(dim(K_gram), c(40L, 40L))
  expect_equal(diag(K_gram), rep(1, 40), tolerance = 1e-6)
  expect_true(all(is.finite(K_gram)))

  K_cos <- hrf_smoothing_kernel(40, TR = 1.5, form = form, method = "cosine",
                                normalise = FALSE, buffer_scans = 2L)
  expect_equal(dim(K_cos), c(40L, 40L))
  expect_true(all(is.finite(K_cos)))

  expect_error(hrf_smoothing_kernel(20, method = "nope"))
})

test_that("meta_fit helpers convert effects and fit contrasts/cov", {
  set.seed(7)
  n <- 12L
  p <- 5L
  Y <- matrix(rnorm(n * p), n, p)
  V <- matrix(runif(n * p, 0.05, 0.2), n, p)
  X <- cbind(1, rnorm(n))

  d1 <- t_to_d(t = 2.5, df = 18)
  expect_equal(names(d1), c("d", "v"))
  expect_true(d1$d > 0)
  d2 <- t_to_d(t = 2.5, df = 18, n = 20)
  expect_true(is.finite(d2$v))

  rz <- r_to_z(0.4, n = 30)
  expect_equal(rz$v, 1 / 27)
  expect_equal(z_to_r(rz$z), 0.4, tolerance = 1e-10)

  n_eff <- meta_effective_n(v = rep(0.05, 3), tau2 = 0.01)
  w <- 1 / (rep(0.05, 3) + 0.01)
  expect_equal(n_eff, sum(w)^2 / sum(w^2))

  Cmat <- matrix(c(0, 1), nrow = 2)
  cres <- fmri_meta_fit_contrasts(Y, V, X, Cmat, method = "fe")
  expect_true(all(c("c_beta", "c_se", "c_z", "ok") %in% names(cres)))
  expect_equal(length(cres$ok), p)

  expect_error(
    fmri_meta_fit_contrasts(Y, V, X, matrix(1, 3, 1), method = "fe"),
    "nrow\\(Cmat\\) must equal ncol\\(X\\)"
  )

  covres <- fmri_meta_fit_cov(Y, V, X, method = "dl")
  expect_true(!is.null(covres$cov_tri))
  expect_equal(ncol(covres$cov_tri), p)

  expect_warning(
    fmri_meta_fit(Y, cbind(V[, 1], -1, V[, 3:p]), X, method = "fe"),
    "Non-positive variances"
  )
})

test_that("fmri_rlm wraps robust fitting and print method", {
  fx <- tiny_matrix_dataset()
  fit <- fmri_rlm(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = fx$dataset,
    strategy = "runwise",
    robust_psi = "huber"
  )
  expect_s3_class(fit, "fmri_rlm")
  expect_s3_class(fit, "fmri_lm")

  out <- capture.output(print(fit))
  expect_true(any(grepl("Robust fMRI Linear Model", out)))

  expect_error(
    fmri_rlm(
      onset ~ hrf(condition),
      block = ~ run,
      dataset = fx$dataset,
      robust_psi = "not-a-psi"
    )
  )
})

test_that("autoplot.Reg and correlation_map cover plotting helpers", {
  skip_if_not_installed("ggplot2")

  reg <- fmrihrf::regressor(
    onsets = c(5, 15, 25),
    hrf = fmrihrf::HRF_SPMG1,
    duration = 0,
    amplitude = 1
  )
  p1 <- ggplot2::autoplot(reg)
  expect_s3_class(p1, "ggplot")

  # Multi-basis HRF path
  reg2 <- fmrihrf::regressor(
    onsets = c(5, 20),
    hrf = fmrihrf::HRF_SPMG3,
    duration = 0,
    amplitude = 1
  )
  p2 <- ggplot2::autoplot(reg2, grid = seq(0, 40, by = 0.5))
  expect_s3_class(p2, "ggplot")

  # Empty-onset / null-grid path
  reg_empty <- fmrihrf::regressor(
    onsets = numeric(0),
    hrf = fmrihrf::HRF_SPMG1,
    duration = 0,
    amplitude = 1
  )
  p3 <- ggplot2::autoplot(reg_empty)
  expect_s3_class(p3, "ggplot")

  fx <- tiny_matrix_dataset()
  bmod <- baseline_model(basis = "poly", degree = 2, sframe = fx$dataset$sampling_frame)
  cm <- correlation_map(bmod, half_matrix = TRUE, method = "pearson")
  expect_s3_class(cm, "ggplot")
  cm2 <- correlation_map(bmod, half_matrix = FALSE, method = "spearman",
                         absolute_limits = FALSE)
  expect_s3_class(cm2, "ggplot")

  expect_error(
    fmrireg:::.correlation_map_common(matrix(1:4, 2, 1)),
    "ncol\\(DM\\) >= 2"
  )
})

test_that("design_plot builds shiny app and validates term_name", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("bslib")
  skip_if_not_installed("thematic")

  fx <- tiny_matrix_dataset()
  emod <- event_model(
    onset ~ hrf(condition),
    data = fx$events,
    block = ~ run,
    sampling_frame = fx$dataset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = fx$dataset$sampling_frame)
  fmod <- fmri_model(emod, bmod, fx$dataset)

  app <- design_plot(fmod, longnames = TRUE, plot_title = "Design")
  expect_s3_class(app, "shiny.appobj")

  expect_error(design_plot(fmod, term_name = "does_not_exist"), "term_name not found")
  expect_error(design_plot(list()), "fmri_model")
})

test_that("mixed_solve_cpp and instantaneous_correlation_rcpp cover thin wrappers", {
  set.seed(1)
  n <- 40
  y <- rnorm(n)
  Z <- cbind(1, rnorm(n))
  X <- matrix(1, n, 1)
  out <- mixed_solve_cpp(y, Z = Z, X = X, SE = TRUE, return_Hinv = TRUE)
  expect_true(is.list(out))
  expect_true(all(c("Vu", "Ve", "beta") %in% names(out)))

  x <- rnorm(50)
  y2 <- 0.5 * x + rnorm(50, sd = 0.2)
  rho <- instantaneous_correlation_rcpp(x, y2, tau_half = 5, fill = "zero")
  expect_equal(length(rho), 50)
  expect_true(any(is.finite(rho)))

  rho_na <- instantaneous_correlation_rcpp(x, y2, eta = 0.2, offset = 2L, fill = "na")
  expect_equal(length(rho_na), 50)

  rho_last <- instantaneous_correlation_rcpp(x, y2, eta = 0.2, offset = -1L, fill = "last")
  expect_equal(length(rho_last), 50)

  expect_error(instantaneous_correlation_rcpp(x, y2[1:10], eta = 0.1))
  expect_error(instantaneous_correlation_rcpp(x, y2))
  expect_error(instantaneous_correlation_rcpp(x, y2, eta = -1))
})

test_that("lss_fast internal helper and deprecated lss_compute_r error", {
  fx <- tiny_matrix_dataset(n_time = 80, n_vox = 4, n_events = 8)
  fit_prep <- estimate_betas(
    fx$dataset,
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "ols"
  )
  # Direct unit exercise of internal LSS helper via a small synthetic bdes
  Y <- as.matrix(fx$dataset$datamat)
  dmat <- as.matrix(design_matrix(
    event_model(
      onset ~ hrf(condition),
      data = fx$events,
      block = ~ run,
      sampling_frame = fx$dataset$sampling_frame
    )
  ))
  bdes <- list(
    dmat_ran = dmat,
    dmat_base = matrix(1, nrow(dmat), 1),
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  betas <- fmrireg:::lss_fast(fx$dataset, bdes, Y = Y, use_cpp = TRUE)
  expect_true(is.matrix(betas))
  expect_equal(nrow(betas), ncol(dmat))

  expect_error(
    fmrireg:::lss_fast(fx$dataset, bdes, Y = Y[1:10, , drop = FALSE]),
    "must have"
  )
  expect_error(fmrireg:::lss_compute_r(NULL, NULL), "deprecated")
})
