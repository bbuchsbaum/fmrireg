# Pass 3, stream R2-2: robust = "t" is rejected everywhere with a clear error
# (it was accepted and silently fitted as "none" or "huber"), and combine = on
# beta + SE input without t-statistics errors clearly instead of failing inside
# the fitter (nifti/h5) or silently running a meta-regression labelled
# "combine:..." (GDS).

mock_nifti_gd_r22 <- function(beta, se, tvals = NULL, df = NULL) {
  S <- nrow(beta)
  P <- ncol(beta)
  dims <- c(P, 1L, 1L)
  store <- new.env()
  bp <- paste0("beta_", seq_len(S))
  sp <- paste0("se_", seq_len(S))
  tp <- if (!is.null(tvals)) paste0("t_", seq_len(S)) else NULL
  for (i in seq_len(S)) {
    store[[bp[i]]] <- array(beta[i, ], dims)
    store[[sp[i]]] <- array(se[i, ], dims)
    if (!is.null(tvals)) store[[tp[i]]] <- array(tvals[i, ], dims)
  }
  gd <- structure(
    list(beta_paths = bp, se_paths = sp, var_paths = NULL, t_paths = tp, df = df,
         subjects = paste0("s", seq_len(S)), covariates = NULL, mask = NULL,
         mask_data = NULL, dim = dims, format = "nifti", n_subjects = S, n_voxels = P),
    class = c("group_data_nifti", "group_data")
  )
  list(gd = gd, store = store)
}

r22_sim <- function(S = 10L, P = 4L) {
  set.seed(22)
  beta <- matrix(stats::rnorm(S * P, 0.3), S, P)
  se <- matrix(stats::runif(S * P, 0.2, 1), S, P)
  list(beta = beta, se = se, S = S, P = P)
}

r22_gds <- function(with_t = FALSE) {
  set.seed(221)
  S <- 8L
  df <- expand.grid(subject = paste0("s", seq_len(S)), roi = c("R1", "R2"),
                    stringsAsFactors = FALSE)
  df$beta <- stats::rnorm(nrow(df), 0.2)
  df$se <- stats::runif(nrow(df), 0.2, 0.6)
  cols <- c(beta = "beta", se = "se")
  if (with_t) {
    df$t <- df$beta / df$se
    df$df <- 30
    cols <- c(cols, t = "t", df = "df")
  }
  group_data(df, format = "csv", effect_cols = cols,
             subject_col = "subject", roi_col = "roi")
}

robust_t_msg <- "robust = \"t\".*not implemented"

test_that("R2-2: fmri_meta rejects robust = 't' (nifti path)", {
  sim <- r22_sim()
  mk <- mock_nifti_gd_r22(sim$beta, sim$se)
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]]
  )
  expect_error(fmri_meta(mk$gd, robust = "t", verbose = FALSE, n_threads = 1),
               robust_t_msg)
  # huber and none still run
  expect_s3_class(fmri_meta(mk$gd, robust = "huber", verbose = FALSE, n_threads = 1),
                  "fmri_meta")
})

test_that("R2-2: fmri_meta rejects robust = 't' (GDS path, with and without contrasts)", {
  skip_if_not_installed("fmrigds")
  gd <- r22_gds()
  expect_error(fmri_meta(gd, robust = "t", verbose = FALSE), robust_t_msg)
  expect_error(fmri_meta(gd, robust = "t", contrasts = c("(Intercept)" = 1),
                         verbose = FALSE), robust_t_msg)
})

test_that("R2-2: fitters, C++ kernels and GDS reducers reject robust = 't'", {
  set.seed(3)
  Y <- matrix(stats::rnorm(20), 10, 2)
  V <- matrix(0.5, 10, 2)
  X <- matrix(1, 10, 1)
  expect_error(fmrireg:::.meta_fit_core(Y, V, X, method = "pm", robust = "t"), robust_t_msg)
  expect_error(meta_fit_cov_cpp(Y, V, X, "pm", "t", 1.345, 2L, 1L), "robust")
  expect_error(meta_fit_cpp(Y, V, X, "pm", "t", 1.345, 2L, 1L), "robust")
  expect_error(fmri_meta_fit(Y, V, X, robust = "t"), robust_t_msg)
  red <- fmrireg:::.fmrigds_make_meta_reducer("pm")
  expect_error(red(beta = Y, var = V, X = X, opts = list(robust = "t")), robust_t_msg)
  reg <- fmrireg:::.fmrigds_make_meta_reg("pm")
  expect_error(reg(beta = Y, var = V, X = X, opts = list(robust = "t")), robust_t_msg)
})

test_that("R2-2: fmri_ttest meta engine rejects option fmrireg.meta.robust = 't'", {
  sim <- r22_sim()
  gd <- structure(list(blocks = list(list(Y = sim$beta, V = sim$se^2,
                                          covars = data.frame(s = seq_len(sim$S))))),
                  class = "group_data")
  withr::local_options(fmrireg.meta.robust = "t")
  expect_error(fmri_ttest(gd, formula = ~ 1, engine = "meta"), robust_t_msg)
})

combine_msg <- "combine.*t-statistic"

test_that("R2-2: combine with beta + SE input (no t) errors clearly (nifti path)", {
  sim <- r22_sim()
  mk <- mock_nifti_gd_r22(sim$beta, sim$se)
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]]
  )
  expect_error(fmri_meta(mk$gd, combine = "stouffer", weights = "equal",
                         verbose = FALSE, n_threads = 1),
               combine_msg)
})

test_that("R2-2: combine with beta + SE + t input combines the t-statistics (nifti path)", {
  sim <- r22_sim()
  tvals <- sim$beta / sim$se
  mk <- mock_nifti_gd_r22(sim$beta, sim$se, tvals = tvals, df = 25)
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]]
  )
  fit <- fmri_meta(mk$gd, combine = "stouffer", weights = "equal",
                   verbose = FALSE, n_threads = 1)
  ref <- combine_t_statistics(tvals, df = 25, method = "stouffer", weights = "equal")
  expect_identical(fit$method, "combine:stouffer")
  expect_equal(as.numeric(fit$coefficients), as.numeric(ref), tolerance = 1e-12)
})

test_that("R2-2: combine with beta + SE input (no t) errors clearly (fit_meta_chunk)", {
  sim <- r22_sim()
  expect_error(
    fmrireg:::fit_meta_chunk(list(beta = sim$beta, se = sim$se), matrix(1, sim$S, 1),
                             method = "combine:stouffer", robust = "none",
                             weights = "equal", weights_custom = NULL, n_threads = 1,
                             combine = "stouffer"),
    combine_msg
  )
})

test_that("R2-2: combine with beta + SE input (no t) errors clearly (GDS path)", {
  skip_if_not_installed("fmrigds")
  gd <- r22_gds(with_t = FALSE)
  expect_error(fmri_meta(gd, combine = "stouffer", weights = "equal", verbose = FALSE),
               combine_msg)
})

test_that("R2-2: combine with beta + SE input (no t) errors clearly (legacy CSV path)", {
  set.seed(5)
  S <- 6L
  df <- data.frame(subject = paste0("s", seq_len(S)), beta = stats::rnorm(S),
                   se = stats::runif(S, 0.2, 0.5))
  gd <- suppressWarnings(group_data_from_csv(df, effect_cols = c(beta = "beta", se = "se"),
                                             subject_col = "subject"))
  withr::local_options(fmrireg.suppress_deprecation = TRUE)
  expect_error(fmri_meta(gd, combine = "stouffer", weights = "equal", verbose = FALSE),
               combine_msg)
})
