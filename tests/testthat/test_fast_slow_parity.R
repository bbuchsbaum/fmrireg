# Certifies that the fast (matrix) and slow (formula/lm) engines agree, so the
# fast path is a safe default. Betas and t/F statistics must match to 1e-6.
#
# NOTE: for AR fits the two engines intentionally differ on p-values only: the
# fast path adjusts the residual df via calculate_effective_df() (whitening has
# a df cost) while the slow path uses the raw residual df. So `prob` is compared
# only for plain OLS, never for AR.

make_ds <- function(seed = 3, nvox = 4, ar_phi = NULL, n = 150, runs = 1, aliased = FALSE) {
  set.seed(seed)
  per <- n / runs
  ev <- do.call(rbind, lapply(seq_len(runs), function(r) {
    data.frame(
      onset = sort(runif(12, 0, per - 10)),
      cond = factor(rep(c("A", "B", "C"), 4)),
      run = r
    )
  }))
  mk <- function() {
    if (is.null(ar_phi)) {
      rnorm(n)
    } else {
      e <- numeric(n)
      e[1] <- rnorm(1)
      for (t in 2:n) e[t] <- ar_phi * e[t - 1] + rnorm(1)
      e
    }
  }
  Y <- sapply(seq_len(nvox), function(v) mk())
  if (aliased) {
    ev$m <- rnorm(nrow(ev))
    ev$m2 <- ev$m
  }
  fmridataset::matrix_dataset(Y, TR = 1, run_length = rep(per, runs), event_table = ev)
}

fit_both <- function(form, dset, ...) {
  list(
    fast = suppressWarnings(fmri_lm(form, block = ~run, dataset = dset,
      strategy = "runwise", use_fast_path = TRUE, progress = FALSE, ...)),
    slow = suppressWarnings(fmri_lm(form, block = ~run, dataset = dset,
      strategy = "runwise", use_fast_path = FALSE, progress = FALSE, ...))
  )
}

# Pull a contrast's estimate/stat/prob vectors directly from the result table
# (robust to the stats()/coef() accessor element semantics).
get_con <- function(fit, name) {
  ct <- fit$result$contrasts
  row <- ct[ct$name == name, , drop = FALSE]
  if (nrow(row) == 0L) stop("contrast not found: ", name)
  d <- row$data[[1]]
  list(
    estimate = as.numeric(unlist(d$estimate)),
    stat = as.numeric(unlist(d$stat)),
    prob = if (!is.null(d$prob)) as.numeric(unlist(d$prob)) else NULL
  )
}

expect_engine_parity <- function(fits, cname, compare_prob = TRUE) {
  expect_equal(as.matrix(coef(fits$fast, type = "betas")),
    as.matrix(coef(fits$slow, type = "betas")), tolerance = 1e-6)
  cf <- get_con(fits$fast, cname)
  cs <- get_con(fits$slow, cname)
  expect_equal(cs$estimate, cf$estimate, tolerance = 1e-6)
  expect_equal(cs$stat, cf$stat, tolerance = 1e-6)
  if (compare_prob) expect_equal(cs$prob, cf$prob, tolerance = 1e-6)
}

test_that("OLS: fast == slow (betas, t-stat, p-value)", {
  dset <- make_ds()
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AmB"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con), dset)
  expect_engine_parity(fits, "AmB", compare_prob = TRUE)
})

test_that("oneway F-contrast: fast == slow (betas, F-stat)", {
  dset <- make_ds()
  con <- contrast_set(oneway_contrast(~cond, name = "cond_F"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con), dset)
  expect_engine_parity(fits, "cond_F", compare_prob = FALSE)
})

test_that("AR(1): fast == slow betas and t-stat (p-values differ by design)", {
  dset <- make_ds(ar_phi = 0.5)
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AmB"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con), dset, cor_struct = "ar1")
  expect_equal(as.numeric(ar_parameters(fits$fast)), as.numeric(ar_parameters(fits$slow)),
    tolerance = 1e-6)
  expect_engine_parity(fits, "AmB", compare_prob = FALSE)
})

test_that("AR(2): fast == slow betas and t-stat", {
  dset <- make_ds(ar_phi = 0.4)
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AmB"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con), dset, cor_struct = "ar2")
  expect_equal(as.numeric(ar_parameters(fits$fast)), as.numeric(ar_parameters(fits$slow)),
    tolerance = 1e-6)
  expect_engine_parity(fits, "AmB", compare_prob = FALSE)
})

test_that("multi-run pooled: fast == slow (betas, t-stat, p-value)", {
  dset <- make_ds(runs = 2, n = 200)
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AmB"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con), dset)
  expect_engine_parity(fits, "AmB", compare_prob = TRUE)
})

test_that("voxelwise AR dispatches regardless of use_fast_path (decoupled gate)", {
  dset <- make_ds(ar_phi = 0.5)
  form <- onset ~ hrf(cond)
  for (fp in c(TRUE, FALSE)) {
    f <- suppressWarnings(fmri_lm(form,
      block = ~run, dataset = dset, strategy = "runwise",
      use_fast_path = fp, cor_struct = "ar1", ar_voxelwise = TRUE, progress = FALSE
    ))
    expect_false(is.null(ar_parameters(f))) # voxelwise AR actually ran
    expect_true(all(is.finite(as.matrix(coef(f, type = "betas")))))
  }
})

test_that("rank-deficient: estimable contrast agrees and stays finite on both engines", {
  dset <- make_ds(aliased = TRUE)
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "AmB"))
  fits <- fit_both(onset ~ hrf(cond, contrasts = con) + hrf(m) + hrf(m2), dset)
  expect_true(all(is.finite(get_con(fits$fast, "AmB")$stat)) &&
    all(is.finite(get_con(fits$slow, "AmB")$stat)))
  expect_engine_parity(fits, "AmB", compare_prob = TRUE)
})
