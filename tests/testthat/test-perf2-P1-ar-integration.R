# Performance pass 2, item P1-5: AR integration overheads.
# - whiten_glm_context() does not compute OLS residuals it does not use
#   (fixed phi).
# - The legacy whitening wrappers honour getOption("fmrireg.ar.parallel_whiten")
#   and never request threads for a single-column response.

.p15_data <- function(n = 60L, V = 4L, seed = 6L) {
  set.seed(seed)
  X <- cbind(1, rnorm(n), cos(seq_len(n) / 4))
  Y <- vapply(seq_len(V), function(v) as.numeric(arima.sim(list(ar = 0.5), n)),
              numeric(n)) + drop(X %*% c(2, 1, -1))
  list(X = X, Y = Y)
}

.p15_ctx <- function(d) {
  fmrireg:::glm_context(X = d$X, Y = d$Y, proj = fmrireg:::.fast_preproject(d$X))
}

test_that("whiten_glm_context with fixed phi equals direct fmriAR whitening", {
  d <- .p15_data()
  out <- fmrireg:::whiten_glm_context(.p15_ctx(d), list(cor_struct = "ar1", phi = 0.4))
  # Pass 3 (R1-4): exact_first defaults to FALSE (the documented default);
  # this path alone used to default it to TRUE.
  plan <- fmriAR::compat$plan_from_phi(phi = list(0.4), theta = list(numeric(0)),
                                       runs = NULL, pooling = "global",
                                       exact_first = FALSE)
  ref <- fmriAR::whiten_apply(plan, X = d$X, Y = d$Y, parallel = FALSE)
  expect_equal(out$X, ref$X, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(out$Y, ref$Y, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(unname(unlist(out$phi_hat)), 0.4)
})

test_that("whiten_glm_context with estimated phi uses the OLS residuals", {
  d <- .p15_data()
  opts <- list(cor_struct = "ar1")
  out <- fmrireg:::whiten_glm_context(.p15_ctx(d), opts)
  resid <- d$Y - d$X %*% qr.solve(d$X, d$Y)
  plan <- fmrireg:::.estimate_ar_via_fmriAR(
    resid, fmrireg:::.normalize_ar_options(opts), NULL, design = d$X)
  ref <- fmriAR::whiten_apply(plan, X = d$X, Y = d$Y, parallel = FALSE)
  expect_equal(out$X, ref$X, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(out$Y, ref$Y, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("whiten_glm_context with fixed phi forms no OLS fit", {
  d <- .p15_data()
  pinv_calls <- 0L
  orig <- fmrireg:::.preproject_pinv
  testthat::local_mocked_bindings(
    .preproject_pinv = function(...) { pinv_calls <<- pinv_calls + 1L; orig(...) }
  )
  fmrireg:::whiten_glm_context(.p15_ctx(d), list(cor_struct = "ar1", phi = 0.4))
  expect_equal(pinv_calls, 0L)
})

test_that("legacy whitening results do not depend on the threading choice", {
  d <- .p15_data()
  for (fun in c("transform", "inplace")) {
    run <- function() {
      if (fun == "transform") {
        fmrireg:::ar_whiten_transform(d$X, d$Y, c(0.3, 0.1), exact_first = TRUE)
      } else {
        fmrireg:::ar_whiten_inplace(d$Y, d$X, c(0.3, 0.1), exact_first_ar1 = TRUE)
      }
    }
    a <- withr::with_options(list(fmrireg.ar.parallel_whiten = TRUE), run())
    b <- withr::with_options(list(fmrireg.ar.parallel_whiten = FALSE), run())
    expect_equal(a$X, b$X, tolerance = 1e-14)
    expect_equal(a$Y, b$Y, tolerance = 1e-14)
    plan <- fmriAR::compat$plan_from_phi(phi = c(0.3, 0.1),
                                         theta = list(numeric(0)),
                                         pooling = "global", exact_first = TRUE)
    ref <- fmriAR::whiten_apply(plan, X = d$X, Y = d$Y, parallel = FALSE)
    expect_equal(b$X, ref$X, tolerance = 1e-14, ignore_attr = TRUE)
    expect_equal(b$Y, ref$Y, tolerance = 1e-14, ignore_attr = TRUE)
  }
})

test_that("legacy whitening honours fmrireg.ar.parallel_whiten and single columns", {
  d <- .p15_data()
  seen <- logical(0)
  orig <- fmriAR::whiten_apply
  testthat::local_mocked_bindings(
    whiten_apply = function(plan, X, Y, runs = NULL, run_starts = NULL,
                            parallel = TRUE, ...) {
      seen <<- c(seen, parallel)
      orig(plan, X = X, Y = Y, runs = runs, parallel = FALSE, ...)
    },
    .package = "fmriAR"
  )
  calls <- list(
    transform = function(Y) fmrireg:::ar_whiten_transform(d$X, Y, 0.4),
    inplace = function(Y) fmrireg:::ar_whiten_inplace(Y, d$X, 0.4)
  )
  for (nm in names(calls)) {
    f <- calls[[nm]]
    seen <- logical(0)
    withr::with_options(list(fmrireg.ar.parallel_whiten = NULL), f(d$Y))
    expect_identical(seen, TRUE, info = paste(nm, "default"))
    seen <- logical(0)
    withr::with_options(list(fmrireg.ar.parallel_whiten = FALSE), f(d$Y))
    expect_identical(seen, FALSE, info = paste(nm, "option FALSE"))
    seen <- logical(0)
    withr::with_options(list(fmrireg.ar.parallel_whiten = TRUE),
                        f(d$Y[, 1, drop = FALSE]))
    expect_identical(seen, FALSE, info = paste(nm, "single column"))
  }
})
