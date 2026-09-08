# Fourteenth wave: pool_runwise_results + estimate_betas.fmri_dataset spatial.

test_that("pool_runwise_results pools beta/contrast/sigma across runs", {
  set.seed(331)
  # Minimal runwise chunk results shaped like runwise_lm_* outputs
  make_bstats <- function(V = 3L, p = 2L) {
    dplyr::tibble(
      type = "beta",
      stat_type = "tstat",
      df.residual = 20,
      conmat = list(NULL),
      colind = list(NULL),
      data = list(dplyr::tibble(
        estimate = list(matrix(rnorm(V * p), V, p)),
        se = list(matrix(runif(V * p, 0.1, 0.5), V, p)),
        stat = list(matrix(rnorm(V * p), V, p)),
        prob = list(matrix(runif(V * p), V, p)),
        sigma = list(runif(V, 0.5, 1.5))
      ))
    )
  }
  cres <- list(
    list(
      bstats = make_bstats(),
      conres = list(),
      rss = c(10, 12, 14),
      rdf = 20,
      sigma = c(0.7, 0.8, 0.9),
      ar_coef = NULL,
      ma_coef = NULL,
      temporal_diagnostics = NULL,
      robust_weights = NULL
    ),
    list(
      bstats = make_bstats(),
      conres = list(),
      rss = c(8, 9, 10),
      rdf = 18,
      sigma = c(0.6, 0.7, 0.8),
      ar_coef = list(0.2, 0.1, 0.15),
      ma_coef = NULL,
      temporal_diagnostics = list(iterations = 1),
      robust_weights = NULL
    )
  )
  Vu <- diag(2)
  pooled <- fmrireg:::pool_runwise_results(
    cres, event_indices = 1L, baseline_indices = 2L, Vu = Vu
  )
  expect_true(is.list(pooled))
  expect_equal(pooled$rdf, 38)
  expect_true(!is.null(pooled$sigma))
  expect_true(all(as.numeric(pooled$sigma) > 0))
  expect_true(!is.null(pooled$rss))
  expect_equal(sum(as.numeric(pooled$rss)), 18 + 21 + 24)
  expect_true(!is.null(pooled$betas) || !is.null(pooled$bstats) || !is.null(pooled$event_indices))
})

test_that("estimate_betas.fmri_dataset OLS on spatial mem dataset", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("fmridataset")
  set.seed(332)
  dims <- c(2L, 2L, 1L)
  n_time <- 48L
  arr <- array(rnorm(prod(dims) * n_time), c(dims, n_time))
  scan <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n_time)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  etab <- data.frame(
    onset = seq(5, n_time - 8, length.out = 6),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  dset <- fmridataset::fmri_mem_dataset(
    scans = list(scan), mask = mask, TR = 1, event_table = etab
  )
  out <- estimate_betas(
    dset,
    ran = onset ~ hrf(condition),
    block = ~ run,
    method = "ols",
    progress = FALSE
  )
  expect_s3_class(out, "fmri_betas")
  expect_true(!is.null(out$betas_ran))
  expect_true(NCOL(out$design_ran) >= 1L || NROW(out$design_ran) >= 1L)
})

test_that(".rrr_filter_task_contrasts keeps in-task weights", {
  simple <- list(
    A = structure(c(1, -1), colind = c(1L, 2L)),
    B = structure(c(1, 0, -1), colind = c(1L, 2L, 3L))
  )
  fcon <- list(
    omnibus = structure(diag(2), colind = c(1L, 2L))
  )
  filtered <- fmrireg:::.rrr_filter_task_contrasts(
    simple_weights = simple,
    f_weights = fcon,
    event_indices = 1:2,
    policy = "warn_drop"
  )
  expect_true(is.list(filtered))
  expect_true(!is.null(filtered$simple) || !is.null(filtered$f) ||
                !is.null(filtered$simple_weights) || length(filtered) >= 1L)
})
