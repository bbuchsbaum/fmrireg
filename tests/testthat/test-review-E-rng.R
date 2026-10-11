# Review 2026-10, stream E3: seeded functions must not clobber the global RNG.

.expect_rng_untouched <- function(expr) {
  set.seed(2024)
  before <- get(".Random.seed", envir = globalenv())
  force(expr)
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
}

test_that("E3: simulate_noise_vector(seed=) leaves .Random.seed unchanged", {
  a <- NULL
  .expect_rng_untouched(a <- simulate_noise_vector(50, seed = 7))
  # still reproducible
  expect_equal(a, simulate_noise_vector(50, seed = 7))
})

test_that("E3: simulate_simple_dataset(seed=) leaves .Random.seed unchanged", {
  a <- NULL
  .expect_rng_untouched(a <- simulate_simple_dataset(ncond = 2, nreps = 3, seed = 7))
  expect_equal(a$noisy, simulate_simple_dataset(ncond = 2, nreps = 3, seed = 7)$noisy)
})

test_that("E3: simulate_fmri_matrix(random_seed=) leaves .Random.seed unchanged", {
  a <- NULL
  .expect_rng_untouched(
    a <- simulate_fmri_matrix(n = 2, total_time = 60, TR = 2, noise_type = "white",
                              random_seed = 7)
  )
  b <- simulate_fmri_matrix(n = 2, total_time = 60, TR = 2, noise_type = "white",
                            random_seed = 7)
  expect_equal(a$time_series$datamat, b$time_series$datamat)
})

test_that("E3: .rng_snapshot restores a missing .Random.seed to missing", {
  genv <- globalenv()
  set.seed(1)
  saved <- get(".Random.seed", envir = genv)
  on.exit(assign(".Random.seed", saved, envir = genv), add = TRUE)
  rm(".Random.seed", envir = genv)
  restore <- fmrireg:::.rng_snapshot()
  set.seed(3)
  restore()
  expect_false(exists(".Random.seed", envir = genv, inherits = FALSE))
})
