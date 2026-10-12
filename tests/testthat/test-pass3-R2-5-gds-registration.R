# Pass 3, stream R2-5: fmrireg's fmrigds reducers return exactly the assay
# fields they declare in `provides` (the packed coefficient covariance of the
# regression reducers travels as fmrigds' cov_tri attachment, always); the
# bindings are registered once per session; registration failures surface
# as warnings instead of being swallowed.

skip_if_not_installed("fmrigds")

r25_reset <- function() {
  st <- get0(".fmrigds_reg_state", envir = asNamespace("fmrireg"), inherits = FALSE)
  if (!is.null(st)) st$done <- FALSE
}

r25_inputs <- function() {
  set.seed(25)
  S <- 9L
  P <- 4L
  X <- cbind("(Intercept)" = 1, age = stats::rnorm(S))
  beta <- matrix(stats::rnorm(S * P), S, P) + 0.3
  var <- matrix(stats::runif(S * P, 0.05, 0.4), S, P)
  list(S = S, P = P, X = X, beta = beta, var = var)
}

test_that("R2-5: every registered reducer returns exactly its declared fields, every branch", {
  r25_reset()
  fmrireg:::.ensure_fmrigds_registered()
  d <- r25_inputs()
  names_all <- c("meta:fe", "meta:re:dl", "meta:re:pm", "meta:re:reml",
                 "frg:meta:fe", "frg:meta:re:dl", "frg:meta:re:pm", "frg:meta:re:reml",
                 "meta:fe_reg", "meta:re_reg")
  for (nm in names_all) {
    r <- fmrigds::get_reducer(nm)
    expect_false(is.null(r), info = nm)
    is_reg <- grepl("_reg$", nm)
    expect_false("cov_tri" %in% r$provides, info = nm)  # an attachment, not an assay
    for (method in list(NULL, "fe", "dl", "pm")) {
      for (rc in list(NULL, "tri")) {
        opts <- list(method = method, return_cov = rc)
        X <- if (is_reg) d$X else d$X[, 1, drop = FALSE]
        out <- suppressMessages(r$fun(d$beta, d$var, X, NULL, NULL, NULL, NULL, NULL, opts))
        info <- paste(nm, method %||% "default", rc %||% "nocov")
        expect_setequal(setdiff(names(out), "cov_tri"), r$provides)
        for (f in r$provides) {
          expect_false(is.null(out[[f]]), info = paste(info, f))
        }
        if (is_reg) {
          expect_true(is.matrix(out$cov_tri), info = info)
          expect_equal(dim(out$cov_tri), c(3L, d$P), info = info)
        } else {
          expect_null(out$cov_tri, info = info)
        }
      }
    }
    # fit-time contrasts cannot be carried by fmrigds reducer outputs
    expect_error(
      r$fun(d$beta, d$var, d$X, NULL, NULL, NULL, NULL, NULL,
            list(contrasts = matrix(c(0, 1), ncol = 1))),
      "contrasts", info = nm
    )
  }
})

test_that("R2-5: meta:re_reg values and covariance match the direct fit", {
  r25_reset()
  fmrireg:::.ensure_fmrigds_registered()
  d <- r25_inputs()
  r <- fmrigds::get_reducer("meta:re_reg")
  out <- r$fun(d$beta, d$var, d$X, NULL, NULL, NULL, NULL, NULL, list(method = "dl"))
  ref <- fmrireg:::.meta_fit_core(d$beta, d$var, d$X, method = "dl")
  expect_equal(out$coef, ref$beta, tolerance = 1e-12)
  expect_equal(out$tau2, as.numeric(ref$tau2), tolerance = 1e-12)
  expect_equal(out$cov_tri, ref$cov_tri, tolerance = 1e-12)
})

test_that("R2-5: registration runs once per session across repeated GDS calls", {
  set.seed(251)
  df <- expand.grid(subject = paste0("s", 1:6), roi = c("R1", "R2"),
                    stringsAsFactors = FALSE)
  df$beta <- stats::rnorm(nrow(df), 0.2)
  df$se <- 0.3
  gd <- group_data(df, format = "csv", effect_cols = c(beta = "beta", se = "se"),
                   subject_col = "subject", roi_col = "roi")
  r25_reset()
  n_calls <- 0L
  real <- fmrireg:::.register_fmrireg_bindings
  testthat::local_mocked_bindings(
    .register_fmrireg_bindings = function(...) {
      n_calls <<- n_calls + 1L
      real(...)
    }
  )
  for (i in 1:3) fmri_meta(gd, verbose = FALSE)
  for (i in 1:2) suppressMessages(fmri_ttest(gd, formula = ~ 1))
  expect_equal(n_calls, 1L)
})

test_that("R2-5: registration failures surface as warnings and are retried", {
  r25_reset()
  testthat::local_mocked_bindings(
    register_reducer = function(...) stop("registry is read-only"),
    .package = "fmrigds"
  )
  expect_warning(fmrireg:::.ensure_fmrigds_registered(), "registry is read-only")
  # not marked as done: the next call tries again (and warns again)
  expect_warning(fmrireg:::.ensure_fmrigds_registered(), "registry is read-only")
})
