# options(fmrireg.num_threads) / FMRIREG_NUM_THREADS choose the OpenMP thread
# count of the compiled kernels, read at call time; the option wins.

test_that(".fmrireg_num_threads: option wins over the environment variable", {
  withr::local_envvar(FMRIREG_NUM_THREADS = "3")
  withr::local_options(fmrireg.num_threads = 2)
  expect_identical(fmrireg:::.fmrireg_num_threads(), 2L)
  # an explicit 0 option means the OpenMP default even if the env var is set
  withr::local_options(fmrireg.num_threads = 0)
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
  # unusable option values fall back to 0
  withr::local_options(fmrireg.num_threads = -1)
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
  withr::local_options(fmrireg.num_threads = "many")
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
})

test_that(".fmrireg_num_threads: env var is used when the option is unset", {
  withr::local_options(fmrireg.num_threads = NULL)
  withr::local_envvar(FMRIREG_NUM_THREADS = "3")
  expect_identical(fmrireg:::.fmrireg_num_threads(), 3L)
  withr::local_envvar(FMRIREG_NUM_THREADS = "0")
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
  withr::local_envvar(FMRIREG_NUM_THREADS = "abc")
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
  withr::local_envvar(FMRIREG_NUM_THREADS = NA)
  expect_identical(fmrireg:::.fmrireg_num_threads(), 0L)
})

test_that("the thread setting reaches the meta-analysis kernel", {
  seen <- integer(0)
  real <- fmrireg:::meta_fit_cov_cpp
  testthat::local_mocked_bindings(
    meta_fit_cov_cpp = function(Y, V, X, method, robust, huber_c, robust_iter, n_threads) {
      seen <<- c(seen, as.integer(n_threads))
      real(Y, V, X, method, robust, huber_c, robust_iter, n_threads)
    }
  )
  set.seed(1)
  Y <- matrix(rnorm(20), 5, 4); V <- matrix(runif(20, 0.5, 1), 5, 4); X <- matrix(1, 5, 1)

  withr::local_options(fmrireg.num_threads = NULL)
  withr::local_envvar(FMRIREG_NUM_THREADS = "3")
  fit_env <- fmri_meta_fit_cov(Y, V, X, method = "fe")
  expect_identical(seen, 3L)

  withr::local_options(fmrireg.num_threads = 2)
  fit_opt <- fmri_meta_fit_cov(Y, V, X, method = "fe")
  expect_identical(seen, c(3L, 2L))
  invisible(fmrireg:::.meta_fit_core(Y, V, X, method = "fe"))
  expect_identical(seen, c(3L, 2L, 2L))
  # the thread count does not change the result
  expect_identical(fit_opt$beta, fit_env$beta)
  w <- 1 / V
  expect_equal(as.numeric(fit_opt$beta), colSums(w * Y) / colSums(w), tolerance = 1e-12)
})

test_that("loading fmrireg no longer sets RcppParallel's global thread pool", {
  expect_false("RcppParallel" %in% names(getNamespaceImports("fmrireg")))
  expect_false(any(grepl("setThreadOptions", deparse(fmrireg:::.onLoad), fixed = TRUE)))
})
