context("Thread configuration")

test_that(".onLoad sets thread count from option", {
  testthat::local_edition(3)
  skip_if_not_installed("RcppParallel")
  old_threads <- RcppParallel::defaultNumThreads()
  old_opt <- getOption("fmrireg.num_threads")
  on.exit({
    options(fmrireg.num_threads = old_opt)
    RcppParallel::setThreadOptions(numThreads = old_threads)
  })

  options(fmrireg.num_threads = 2)
  fmrireg::.onLoad(NULL, NULL)
  expect_equal(RcppParallel::defaultNumThreads(), 2)
})

 test_that(".onLoad sets thread count from env var", {
   skip_if_not_installed("RcppParallel")
   old_threads <- RcppParallel::defaultNumThreads()
   old_env <- Sys.getenv("FMRIREG_NUM_THREADS", unset = NA)
   on.exit({
     if (is.na(old_env)) Sys.unsetenv("FMRIREG_NUM_THREADS") else Sys.setenv(FMRIREG_NUM_THREADS = old_env)
     RcppParallel::setThreadOptions(numThreads = old_threads)
   })

   Sys.setenv(FMRIREG_NUM_THREADS = 3)
   options(fmrireg.num_threads = NULL)
   fmrireg::.onLoad(NULL, NULL)
   expect_equal(RcppParallel::defaultNumThreads(), 3)
 })
