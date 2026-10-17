# Regression tests for the 2026-10 simplification, stream D (group analysis).

.d_demo_gds <- function() {
  gd_csv <- suppressWarnings(fmrireg:::.demo_group_data_csv())
  fmrireg::group_data(
    gd_csv$data, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )
}

.d_registry_snapshot <- function() {
  red <- c("meta:fe", "meta:re", "meta:fe_reg", "meta:re_reg",
           "meta:re:dl", "meta:re:pm", "meta:re:reml",
           "frg:meta:fe", "frg:meta:re:dl", "frg:meta:re:pm", "frg:meta:re:reml")
  ph <- c("fdr:bh", "fdr:by", "fdr:spatial")
  list(
    reducers = fmrigds::list_reducers(),
    posthoc = fmrigds::list_posthoc(),
    reducer_fun = lapply(stats::setNames(red, red), function(n) {
      r <- fmrigds::get_reducer(n)
      if (is.null(r)) NULL else list(provides = r$provides, fun = r$fun)
    }),
    posthoc_fun = lapply(stats::setNames(ph, ph), function(n) {
      fmrigds::get_posthoc(n)$fun
    })
  )
}

test_that("gds fmri_meta()/fmri_ttest() leave fmrigds's reducer registry alone (B1)", {
  skip_if_not_installed("fmrigds")
  gd <- .d_demo_gds()
  before <- .d_registry_snapshot()
  # fmrigds's own meta:fe provides the *_g outputs
  expect_setequal(before$reducer_fun[["meta:fe"]]$provides,
                  c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"))
  fit <- fmri_meta(gd, formula = ~ 1, method = "pm", verbose = FALSE)
  tt <- fmri_ttest(gd, formula = ~ 1, engine = "meta")
  expect_s3_class(fit, "fmri_meta")
  after <- .d_registry_snapshot()
  expect_identical(after$reducers, before$reducers)
  expect_identical(after$posthoc, before$posthoc)
  for (nm in names(before$reducer_fun)) {
    expect_identical(after$reducer_fun[[nm]], before$reducer_fun[[nm]], info = nm)
  }
  for (nm in names(before$posthoc_fun)) {
    expect_identical(after$posthoc_fun[[nm]], before$posthoc_fun[[nm]], info = nm)
  }
  expect_setequal(fmrigds::get_reducer("meta:fe")$provides,
                  c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"))
  expect_false(any(grepl("^frg:", fmrigds::list_reducers())))
})

test_that("fmrigds meta:fe still computes the closed-form FE estimate after fmrireg runs (B1)", {
  skip_if_not_installed("fmrigds")
  gd <- .d_demo_gds()
  invisible(fmri_meta(gd, formula = ~ 1, method = "fe", verbose = FALSE))
  set.seed(3)
  P <- 4; S <- 6
  b <- array(rnorm(P * S), c(P, S, 1)); v <- array(runif(P * S, 0.5, 2), c(P, S, 1))
  g <- fmrigds::as_gds(list(beta = b, var = v), space = NULL)
  out <- fmrigds::compute(fmrigds::reduce(fmrigds::as_plan(g), method = "meta:fe"))
  est <- as.numeric(fmrigds::assay(out, "beta_g"))
  w <- 1 / v[, , 1]
  expect_equal(est, rowSums(w * b[, , 1]) / rowSums(w), tolerance = 1e-10)
})

test_that("fmrireg's spatial-FDR posthoc is available under a namespaced name", {
  skip_if_not_installed("fmrigds")
  set.seed(42)
  P <- 50; K <- 2
  z_mat <- matrix(rnorm(P * K), nrow = P, ncol = K)
  grp <- rep(1:5, each = 10)
  g <- fmrigds::as_gds(list(z = array(z_mat, dim = c(P, 1, K))), space = NULL)
  pl <- fmrigds::posthoc(fmrigds::as_plan(g), method = "fmrireg:fdr:spatial",
                         options = list(group = grp))
  q <- fmrigds::assay(fmrigds::compute(pl), "q")
  for (k in 1:K) {
    expect_equal(as.numeric(q[, 1, k]), as.numeric(spatial_fdr(z = z_mat[, k], group = grp)$q))
  }
})

test_that("group_data(format = 'fmrilm') fails up front and names the working route (B4)", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  # Previously: "invalid 'file' argument" from the tabular adapter.
  err <- tryCatch(group_data(list(fit, fit), format = "fmrilm"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "group_data(format = \"fmrilm\") is not supported", fixed = TRUE)
  expect_match(conditionMessage(err), "write_results()", fixed = TRUE)
  expect_match(conditionMessage(err), "format = \"nifti\"", fixed = TRUE)
  # A list of fits passed with format = "auto" gets the same message.
  err2 <- tryCatch(group_data(list(fit, fit)), error = function(e) e)
  expect_identical(conditionMessage(err2), conditionMessage(err))
})

# --- B5: contrast weights never silently become zero -------------------------

.d_toy_meta <- function() {
  structure(
    list(
      coefficients = matrix(c(0.3, 0.1, 0.5, -0.2, 0.4, 0.0), nrow = 2,
                            dimnames = list(NULL, c("(Intercept)", "groupB", "age"))),
      se = matrix(c(0.05, 0.06, 0.07, 0.04, 0.05, 0.03), nrow = 2),
      robust = "none"
    ),
    class = "fmri_meta"
  )
}

test_that("contrast.fmri_meta formula: unknown names and unparseable formulas error (B5)", {
  m <- .d_toy_meta()
  # previously all-zero weights (estimate 0, se 0, z NaN) with no message
  expect_error(contrast(m, ~ groupC - age), "groupC")
  expect_error(contrast(m, ~ groupB - groupC), "groupC")
  expect_error(contrast(m, ~ log(groupB)), "Cannot parse contrast formula")
  expect_error(contrast(m, ~ groupB * age), "Cannot parse contrast formula")
})

test_that("contrast.fmri_meta formula weights match the equivalent numeric contrast (B5)", {
  m <- .d_toy_meta()
  ref <- suppressWarnings(contrast(m, c(0, 1, -1)))
  got <- suppressWarnings(contrast(m, ~ groupB - age))
  expect_identical(unname(got$weights), c(0, 1, -1))
  expect_equal(got$estimate, ref$estimate)
  expect_equal(got$se, ref$se)
  expect_identical(unname(contrast(m, ~ groupB)$weights), c(0, 1, 0))
  expect_identical(unname(suppressWarnings(contrast(m, ~ `(Intercept)` + 0.5 * groupB - 2 * age))$weights),
                   c(1, 0.5, -2))
  expect_identical(unname(suppressWarnings(contrast(m, ~ -groupB + age))$weights), c(0, -1, 1))
  expect_identical(unname(suppressWarnings(contrast(m, c(groupB = 1, age = -1)))$weights),
                   c(0, 1, -1))
  expect_error(contrast(m, c(groupB = 1, nope = -1)), "nope")
})

test_that("legacy fit-time contrasts reject unknown coefficient names (B5)", {
  set.seed(5)
  S <- 10L; P <- 4L
  covars <- data.frame(x = rnorm(S))
  beta <- matrix(rnorm(S * P), S, P) + 0.3 * covars$x
  se <- matrix(runif(S * P, 0.2, 1), S, P)
  dims <- c(P, 1L, 1L)
  store <- new.env()
  bp <- paste0("beta_", seq_len(S)); sp <- paste0("se_", seq_len(S))
  for (i in seq_len(S)) {
    store[[bp[i]]] <- array(beta[i, ], dims); store[[sp[i]]] <- array(se[i, ], dims)
  }
  gd <- structure(
    list(beta_paths = bp, se_paths = sp, var_paths = NULL, t_paths = NULL, df = NULL,
         subjects = paste0("s", seq_len(S)), covariates = covars, mask = NULL,
         mask_data = NULL, dim = dims, format = "nifti", n_subjects = S, n_voxels = P),
    class = c("group_data_nifti", "group_data")
  )
  testthat::local_mocked_bindings(read_nifti_data = function(path, mask = NULL) store[[path]])
  fit_named <- fmri_meta(gd, formula = ~ x, method = "fe", verbose = FALSE, n_threads = 1,
                         contrasts = c(x = 1))
  fit_unnamed <- fmri_meta(gd, formula = ~ x, method = "fe", verbose = FALSE, n_threads = 1,
                           contrasts = c(0, 1))
  expect_identical(fit_named$contrasts$estimate, fit_unnamed$contrasts$estimate)
  expect_identical(fit_named$contrasts$se, fit_unnamed$contrasts$se)
  expect_equal(as.numeric(fit_named$contrasts$estimate), as.numeric(fit_named$coefficients[, "x"]))
  # previously the unknown name was appended as a third weight and silently dropped
  expect_error(
    fmri_meta(gd, formula = ~ x, method = "fe", verbose = FALSE, n_threads = 1,
              contrasts = c(z = 1)),
    "z"
  )
  # the h5 method uses the same helper
  expect_identical(fmrireg:::.meta_named_contrast_weights(c(x = 2), c("(Intercept)", "x")),
                   c("(Intercept)" = 0, x = 2))
  expect_error(fmrireg:::.meta_named_contrast_weights(c(z = 1), c("(Intercept)", "x")), "z")
})
