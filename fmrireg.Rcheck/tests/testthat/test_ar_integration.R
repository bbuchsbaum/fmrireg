context("AR modeling integration")

simulate_ar_dataset <- function(ar_coeff = numeric(), n_runs = 2, n_time = 30, n_vox = 4) {
  dat_list <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    run_mat <- replicate(n_vox, {
      if (length(ar_coeff) == 0) {
        rnorm(n_time)
      } else {
        as.numeric(arima.sim(model = list(ar = ar_coeff), n = n_time))
      }
    })
    dat_list[[r]] <- run_mat
  }
  datamat <- do.call(rbind, dat_list)
  event_tab <- expand.grid(run = seq_len(n_runs), onset = c(5, 15))
  event_tab$cond <- factor("A")
  matrix_dataset(datamat, TR = 1, run_length = rep(n_time, n_runs), event_table = event_tab)
}


# Test that AR1 with no correlation matches IID

test_that("iid and ar1 give similar results on white noise", {
  set.seed(1)
  dset <- simulate_ar_dataset(n_runs = 2, ar_coeff = numeric())
  mod_iid <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "iid")
  mod_ar1 <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "ar1")
  expect_equal(coef(mod_iid), coef(mod_ar1), tolerance = 1e-6)

  model <- create_fmri_model(onset ~ hrf(cond), block = ~ run, dataset = dset, durations = 0)
  X <- design_matrix(model)
  proj <- .fast_preproject(X)
  Y <- get_data_matrix(dset)
  ols <- .fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  resid_ols <- Y - ols$fitted
  phi_hat <- .estimate_ar(rowMeans(resid_ols), 1)
  expect_equal(as.numeric(phi_hat), 0, tolerance = 0.1)
})

# Test AR1 recovery and SE comparison

test_that("ar1 recovers phi and adjusts standard errors", {
  set.seed(2)
  phi <- 0.4
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 2)

  model <- create_fmri_model(onset ~ hrf(cond), block = ~ run, dataset = dset, durations = 0)
  X <- design_matrix(model)
  proj <- .fast_preproject(X)
  Y <- get_data_matrix(dset)
  ols <- .fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  resid_ols <- Y - ols$fitted
  phi_hat <- .estimate_ar(rowMeans(resid_ols), 1)
  expect_equal(as.numeric(phi_hat), phi, tolerance = 0.1)

  mod_iid <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "iid")
  mod_ar1 <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "ar1")
  se_iid <- unlist(standard_error(mod_iid))
  se_ar1 <- unlist(standard_error(mod_ar1))
  expect_gt(se_ar1[1], se_iid[1])
})

# Test AR2 recovery

test_that("ar2 recovers coefficients", {
  set.seed(3)
  phi <- c(0.5, -0.25)
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 2)

  model <- create_fmri_model(onset ~ hrf(cond), block = ~ run, dataset = dset, durations = 0)
  X <- design_matrix(model)
  proj <- .fast_preproject(X)
  Y <- get_data_matrix(dset)
  ols <- .fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  resid_ols <- Y - ols$fitted
  phi_hat <- .estimate_ar(rowMeans(resid_ols), 2)
  expect_equal(as.numeric(phi_hat), phi, tolerance = 0.1)
})

# Test global vs runwise

test_that("cor_global gives similar results", {
  set.seed(4)
  phi <- 0.4
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 2)
  mod_run <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "ar1", cor_global = FALSE)
  mod_global <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                        use_fast_path = TRUE, cor_struct = "ar1", cor_global = TRUE)
  expect_equal(coef(mod_run), coef(mod_global), tolerance = 1e-6)
})

# Test ar1_exact_first option

test_that("ar1_exact_first runs", {
  set.seed(5)
  phi <- 0.4
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 1)
  expect_error(
    fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
            use_fast_path = TRUE, cor_struct = "ar1", ar1_exact_first = TRUE),
    NA
  )
})

# Test multiple iterations

test_that("cor_iter > 1 runs", {
  set.seed(6)
  phi <- 0.4
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 1)
  expect_error(
    fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
            use_fast_path = TRUE, cor_struct = "ar1", cor_iter = 2),
    NA
  )
})

# Test ARp recovery for p > 2

test_that("arp recovers coefficients", {
  set.seed(7)
  phi <- c(0.6, -0.3, 0.2)
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 2)

  model <- create_fmri_model(onset ~ hrf(cond), block = ~ run, dataset = dset, durations = 0)
  X <- design_matrix(model)
  proj <- .fast_preproject(X)
  Y <- get_data_matrix(dset)
  ols <- .fast_lm_matrix(X, Y, proj, return_fitted = TRUE)
  resid_ols <- Y - ols$fitted
  phi_hat <- .estimate_ar(rowMeans(resid_ols), length(phi))
  expect_equal(as.numeric(phi_hat), phi, tolerance = 0.1)

  mod_arp <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "arp", ar_p = length(phi))
  expect_true(!is.null(coef(mod_arp)))
})

# Test arp with p=1 matches ar1

test_that("arp with p=1 matches ar1", {
  set.seed(8)
  phi <- 0.5
  dset <- simulate_ar_dataset(ar_coeff = phi, n_runs = 2)
  mod_ar1 <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                     use_fast_path = TRUE, cor_struct = "ar1")
  mod_arp1 <- fmri_lm(onset ~ hrf(cond), block = ~ run, dataset = dset,
                      use_fast_path = TRUE, cor_struct = "arp", ar_p = 1)
  expect_equal(coef(mod_ar1), coef(mod_arp1), tolerance = 1e-6)
})

