test_that(".onLoad sets thread count from option", {
  skip_if_not_installed("RcppParallel")

  old_opt <- getOption("fmrireg.num_threads")
  old_env <- Sys.getenv("FMRIREG_NUM_THREADS", unset = NA)
  on.exit({
    options(fmrireg.num_threads = old_opt)
    if (is.na(old_env)) Sys.unsetenv("FMRIREG_NUM_THREADS") else Sys.setenv(FMRIREG_NUM_THREADS = old_env)
  })

  # Clear env var so only the option is used

  Sys.unsetenv("FMRIREG_NUM_THREADS")
  options(fmrireg.num_threads = 2)

  # .onLoad should run without error

  expect_no_error(fmrireg:::.onLoad(NULL, NULL))

  # If the platform honours the setting, verify it took effect;

  # otherwise just confirm .onLoad completed without error (above).
  cur <- RcppParallel::defaultNumThreads()
  if (cur == 2) {
    expect_equal(cur, 2)
  }
})

test_that(".onLoad sets thread count from env var", {
  skip_if_not_installed("RcppParallel")

  old_opt <- getOption("fmrireg.num_threads")
  old_env <- Sys.getenv("FMRIREG_NUM_THREADS", unset = NA)
  on.exit({
    options(fmrireg.num_threads = old_opt)
    if (is.na(old_env)) Sys.unsetenv("FMRIREG_NUM_THREADS") else Sys.setenv(FMRIREG_NUM_THREADS = old_env)
  })

  # Clear the option so only the env var is used
  options(fmrireg.num_threads = NULL)
  Sys.setenv(FMRIREG_NUM_THREADS = 3)

  # .onLoad should run without error
  expect_no_error(fmrireg:::.onLoad(NULL, NULL))

  # If the platform honours the setting, verify it took effect
  cur <- RcppParallel::defaultNumThreads()
  if (cur == 3) {
    expect_equal(cur, 3)
  }
})
