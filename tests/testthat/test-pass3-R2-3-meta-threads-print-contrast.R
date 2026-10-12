# Pass 3, stream R2-3: fmri_meta forwards the resolved thread count to its
# methods; print.fmri_meta reports I^2 (a fraction) as a percentage; formula
# contrasts error on unknown coefficient names instead of returning all-zero
# weights.

mock_nifti_gd_r23 <- function(beta, se) {
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
         subjects = paste0("s", seq_len(S)), covariates = NULL, mask = NULL,
         mask_data = NULL, dim = dims, format = "nifti", n_subjects = S, n_voxels = P),
    class = c("group_data_nifti", "group_data")
  )
  list(gd = gd, store = store)
}

test_that("R2-3: fmri_meta methods receive the resolved thread count", {
  set.seed(31)
  beta <- matrix(stats::rnorm(24), 8, 3)
  se <- matrix(0.5, 8, 3)
  mk <- mock_nifti_gd_r23(beta, se)
  seen <- integer(0)
  real_chunk <- fmrireg:::fit_meta_chunk
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) mk$store[[path]],
    fit_meta_chunk = function(...) {
      seen <<- c(seen, list(...)$n_threads)
      real_chunk(...)
    }
  )
  expected_auto <- max(1L, parallel::detectCores() - 1L)

  withr::local_options(fmrireg.num_threads = 0)
  fmri_meta(mk$gd, verbose = FALSE)
  expect_equal(seen, expected_auto)

  seen <- integer(0)
  fmri_meta(mk$gd, verbose = FALSE, n_threads = 0)
  expect_equal(seen, expected_auto)

  seen <- integer(0)
  withr::local_options(fmrireg.num_threads = 2)
  fmri_meta(mk$gd, verbose = FALSE)
  expect_equal(seen, 2)

  seen <- integer(0)
  fmri_meta(mk$gd, verbose = FALSE, n_threads = 1)
  expect_equal(seen, 1)
})

test_that("R2-3: GDS fmri_meta receives the resolved thread count", {
  skip_if_not_installed("fmrigds")
  set.seed(32)
  df <- expand.grid(subject = paste0("s", 1:6), roi = c("R1", "R2"),
                    stringsAsFactors = FALSE)
  df$beta <- stats::rnorm(nrow(df))
  df$se <- 0.4
  gd <- group_data(df, format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject", roi_col = "roi")
  seen <- integer(0)
  real <- fmrireg:::.gds_meta_fit_chunked
  testthat::local_mocked_bindings(
    .gds_meta_fit_chunked = function(...) {
      seen <<- c(seen, list(...)$n_threads)
      real(...)
    }
  )
  withr::local_options(fmrireg.num_threads = 0)
  fmri_meta(gd, verbose = FALSE)
  expect_equal(unique(seen), max(1L, parallel::detectCores() - 1L))
})

test_that("R2-3: print.fmri_meta reports I^2 as a percentage of the fraction", {
  fit <- structure(
    list(coefficients = matrix(0.2, nrow = 2), se = matrix(0.05, nrow = 2),
         method = "pm", robust = "none", formula = ~ 1, n_subjects = 10,
         n_voxels = 2, tau2 = c(0.1, 0.3), I2 = c(0.25, 0.5)),
    class = "fmri_meta"
  )
  out <- capture.output(print(fit))
  i2_line <- grep("I\\^2", out, value = TRUE)
  expect_length(i2_line, 1L)
  expect_match(i2_line, "37.5%", fixed = TRUE)
  expect_false(grepl("0.375", i2_line, fixed = TRUE))
})

toy_meta_r23 <- function() {
  structure(
    list(coefficients = matrix(c(0.3, 0.1, 0.2, 0.5, 0.0, -0.1), nrow = 2,
                               dimnames = list(NULL, c("(Intercept)", "groupB", "age"))),
         se = matrix(0.1, nrow = 2, ncol = 3),
         robust = "none"),
    class = "fmri_meta"
  )
}

test_that("R2-3: parse_contrast_formula errors naming unmatched terms", {
  x <- toy_meta_r23()
  expect_error(fmrireg:::parse_contrast_formula(~ groupold - groupyoung, x),
               "groupold.*groupyoung")
  expect_error(fmrireg:::parse_contrast_formula(~ groupB - agee, x), "agee")
  expect_error(contrast(x, ~ groupb), "groupb")
})

test_that("R2-3: parse_contrast_formula builds weights for valid formulas", {
  x <- toy_meta_r23()
  expect_equal(fmrireg:::parse_contrast_formula(~ groupB - age, x),
               c("(Intercept)" = 0, groupB = 1, age = -1))
  expect_equal(fmrireg:::parse_contrast_formula(~ groupB, x),
               c("(Intercept)" = 0, groupB = 1, age = 0))
  expect_equal(fmrireg:::parse_contrast_formula(~ `(Intercept)` + 2 * age, x),
               c("(Intercept)" = 1, groupB = 0, age = 2))
  con <- expect_warning(contrast(x, ~ groupB - age), "approximate")
  expect_equal(con$estimate, as.numeric(x$coefficients %*% c(0, 1, -1)))
})
