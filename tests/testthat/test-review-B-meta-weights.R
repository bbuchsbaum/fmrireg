# Review 2026-10, stream B: weights = "equal" / "custom" are relative weights
# with an estimated residual scale (B2), and custom weight matrices are sliced
# per chunk.

mock_nifti_gd_w <- function(beta, se, covariates = NULL) {
  S <- nrow(beta)
  P <- ncol(beta)
  dims <- c(P, 1L, 1L)
  store <- new.env()
  bp <- paste0("beta_", seq_len(S))
  sp <- paste0("se_", seq_len(S))
  for (i in seq_len(S)) {
    store[[bp[i]]] <- array(beta[i, ], dims)
    store[[sp[i]]] <- array(se[i, ], dims)
  }
  gd <- structure(
    list(beta_paths = bp, se_paths = sp, var_paths = NULL, t_paths = NULL, df = NULL,
         subjects = paste0("s", seq_len(S)), covariates = covariates, mask = NULL,
         mask_data = NULL, dim = dims, format = "nifti", n_subjects = S, n_voxels = P),
    class = c("group_data_nifti", "group_data")
  )
  list(gd = gd, store = store)
}

weights_sim <- function(S = 13L, P = 9L, scale = 1, seed = 21) {
  set.seed(seed)
  covars <- data.frame(x = stats::rnorm(S))
  beta <- scale * (matrix(stats::rnorm(S * P), S, P) + 0.4 + 0.3 * covars$x)
  se <- matrix(stats::runif(S * P, 0.1, 2), S, P)  # deliberately unrelated to the data
  list(beta = beta, se = se, covars = covars)
}

test_that("B2: equal-weights fit equals the one-sample t-test for any data scale", {
  for (scale in c(1e-3, 1, 250)) {
    sim <- weights_sim(scale = scale)
    mk <- mock_nifti_gd_w(sim$beta, sim$se)
    testthat::local_mocked_bindings(
      read_nifti_data = function(path, mask = NULL) mk$store[[path]]
    )
    for (method in c("fe", "pm")) {
      fit <- fmri_meta(mk$gd, method = method, weights = "equal", verbose = FALSE, n_threads = 1)
      for (j in seq_len(ncol(sim$beta))) {
        tt <- stats::t.test(sim$beta[, j])
        info <- paste("scale", scale, method, "voxel", j)
        expect_equal(unname(fit$coefficients[j, 1]), unname(tt$estimate), tolerance = 1e-10, info = info)
        expect_equal(unname(fit$se[j, 1]), unname(tt$stderr), tolerance = 1e-10, info = info)
        expect_equal(unname(fit$coefficients[j, 1] / fit$se[j, 1]), unname(tt$statistic),
                     tolerance = 1e-10, info = info)
      }
      expect_true(all(is.na(fit$tau2)))
    }
  }
})

test_that("B2: custom weights equal weighted least squares (lm with weights)", {
  sim <- weights_sim()
  mk <- mock_nifti_gd_w(sim$beta, sim$se, covariates = sim$covars)
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]]
  )
  S <- nrow(sim$beta)
  w <- seq(0.5, 3, length.out = S)
  fit <- fmri_meta(mk$gd, formula = ~ x, weights = "custom", weights_custom = w,
                   verbose = FALSE, n_threads = 1)
  # rescaling the weights does not change anything
  fit10 <- fmri_meta(mk$gd, formula = ~ x, weights = "custom", weights_custom = 10 * w,
                     verbose = FALSE, n_threads = 1)
  expect_equal(fit10$se, fit$se, tolerance = 1e-10)
  for (j in seq_len(ncol(sim$beta))) {
    m <- stats::lm(sim$beta[, j] ~ sim$covars$x, weights = w)
    cf <- summary(m)$coefficients
    expect_equal(unname(fit$coefficients[j, ]), unname(cf[, "Estimate"]), tolerance = 1e-10)
    expect_equal(unname(fit$se[j, ]), unname(cf[, "Std. Error"]), tolerance = 1e-10)
    # stored covariance equals vcov(lm)
    V <- stats::vcov(m)
    expect_equal(fit$cov$tri[, j], c(V[1, 1], V[1, 2], V[2, 2]), tolerance = 1e-10)
  }
})

test_that("B2: a custom weight matrix is sliced per chunk (chunk_size does not matter)", {
  sim <- weights_sim(P = 11L)
  mk <- mock_nifti_gd_w(sim$beta, sim$se, covariates = sim$covars)
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]]
  )
  S <- nrow(sim$beta)
  P <- ncol(sim$beta)
  set.seed(3)
  W <- matrix(stats::runif(S * P, 0.2, 4), S, P)
  fits <- lapply(c(3L, 4L, 100L), function(cs) {
    fmri_meta(mk$gd, formula = ~ x, weights = "custom", weights_custom = W,
              chunk_size = cs, verbose = FALSE, n_threads = 1)
  })
  expect_equal(fits[[2]]$coefficients, fits[[1]]$coefficients)
  expect_equal(fits[[3]]$coefficients, fits[[1]]$coefficients)
  expect_equal(fits[[3]]$se, fits[[1]]$se)
  for (j in seq_len(P)) {
    m <- stats::lm(sim$beta[, j] ~ sim$covars$x, weights = W[, j])
    cf <- summary(m)$coefficients
    expect_equal(unname(fits[[1]]$coefficients[j, ]), unname(cf[, "Estimate"]), tolerance = 1e-10)
    expect_equal(unname(fits[[1]]$se[j, ]), unname(cf[, "Std. Error"]), tolerance = 1e-10)
  }
})

test_that("B2: fmri_ttest meta engine with equal / custom weights matches t.test / lm", {
  sim <- weights_sim(scale = 40)
  gd <- structure(
    list(blocks = list(list(Y = sim$beta, V = sim$se^2, covars = sim$covars))),
    class = "group_data"
  )
  fit_eq <- fmri_ttest(gd, formula = ~ 1, engine = "meta", weights = "equal")
  for (j in seq_len(ncol(sim$beta))) {
    tt <- stats::t.test(sim$beta[, j])
    expect_equal(unname(fit_eq$beta[1, j]), unname(tt$estimate), tolerance = 1e-10)
    expect_equal(unname(fit_eq$se[1, j]), unname(tt$stderr), tolerance = 1e-10)
  }
  w <- seq(1, 2, length.out = nrow(sim$beta))
  fit_cw <- fmri_ttest(gd, formula = ~ 1 + x, engine = "meta", weights = "custom",
                       weights_custom = w, contrast = c(`(Intercept)` = 1, x = 1))
  for (j in seq_len(ncol(sim$beta))) {
    m <- stats::lm(sim$beta[, j] ~ sim$covars$x, weights = w)
    cf <- summary(m)$coefficients
    expect_equal(unname(fit_cw$se[, j]), unname(cf[, "Std. Error"]), tolerance = 1e-10)
    L <- c(1, 1)
    expect_equal(fit_cw$se_contrast[j], sqrt(drop(t(L) %*% stats::vcov(m) %*% L)),
                 tolerance = 1e-10)
  }
})

test_that("B2: fit_meta_single (ROI path) uses the same weighting", {
  sim <- weights_sim()
  X <- stats::model.matrix(~ x, sim$covars)
  w <- seq(0.5, 3, length.out = nrow(X))
  single <- fmrireg:::fit_meta_single(sim$beta[, 1], sim$se[, 1], X, method = "pm",
                                      robust = "none", weights = "custom", weights_custom = w)
  m <- stats::lm(sim$beta[, 1] ~ sim$covars$x, weights = w)
  expect_equal(single$coefficients, unname(stats::coef(m)), tolerance = 1e-10)
  expect_equal(single$se, unname(sqrt(diag(stats::vcov(m)))), tolerance = 1e-10)
})
