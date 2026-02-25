## Parity matrix for rrr_gls:
## - noise model: iid / ar1
## - run structure: single / multi-run
## - rank mode: full rank parity vs reduced-rank approximation
## plus weak-parity sanity against latent_sketch.

.ensure_rrr_engine_registered <- function() {
  if (exists(".register_builtin_engines",
             envir = asNamespace("fmrireg"),
             mode = "function",
             inherits = FALSE)) {
    fmrireg:::.register_builtin_engines()
  }
}

.simulate_rrr_parity_dataset <- function(seed = 1L,
                                         nruns = 1L,
                                         ar1 = FALSE,
                                         T_run = 60L,
                                         V = 48L,
                                         signal_sd = 0.35,
                                         noise_sd = 0.9,
                                         rho = 0.35) {
  set.seed(seed)

  nruns <- as.integer(nruns)
  run_length <- rep(as.integer(T_run), nruns)
  Ttot <- sum(run_length)

  per_run_events <- data.frame(
    onsets = c(6, 18, 30, 42),
    condition = factor(c("A", "B", "A", "B"), levels = c("A", "B"))
  )
  ev <- do.call(
    rbind,
    lapply(seq_len(nruns), function(r) {
      data.frame(onsets = per_run_events$onsets,
                 condition = per_run_events$condition,
                 run = r)
    })
  )

  Y_tmp <- matrix(rnorm(Ttot * V, sd = 0.1), nrow = Ttot, ncol = V)
  dtmp <- fmridataset::matrix_dataset(
    Y_tmp,
    TR = 2,
    run_length = run_length,
    event_table = ev
  )
  model <- create_fmri_model(
    formula = onsets ~ hrf(condition),
    block = ~run,
    dataset = dtmp
  )
  X <- as.matrix(design_matrix(model))

  tm <- fmridesign::term_matrices(model)
  event_indices <- as.integer(attr(tm, "event_term_indices"))

  B_true <- matrix(0, nrow = ncol(X), ncol = V)
  B_true[event_indices, ] <- matrix(
    rnorm(length(event_indices) * V, sd = signal_sd),
    nrow = length(event_indices)
  )
  signal <- X %*% B_true

  noise <- matrix(0, nrow = Ttot, ncol = V)
  if (!isTRUE(ar1)) {
    noise[] <- rnorm(Ttot * V, sd = noise_sd)
  } else {
    row_start <- 1L
    for (r in seq_len(nruns)) {
      n_r <- run_length[r]
      idx <- seq.int(row_start, row_start + n_r - 1L)
      for (v in seq_len(V)) {
        e <- numeric(n_r)
        e[1] <- rnorm(1, sd = noise_sd / sqrt(1 - rho^2))
        for (t in 2:n_r) {
          e[t] <- rho * e[t - 1L] + rnorm(1, sd = noise_sd)
        }
        noise[idx, v] <- e
      }
      row_start <- row_start + n_r
    }
  }

  Y <- signal + noise
  dset <- fmridataset::matrix_dataset(
    Y,
    TR = 2,
    run_length = run_length,
    event_table = ev
  )

  list(dataset = dset, events = ev)
}

.event_beta_matrix <- function(fit) {
  t(as.matrix(fit$result$betas$data[[1]]$estimate[[1]]))
}

test_that("rrr_gls parity matrix against standard fmri_lm", {
  skip_on_cran()
  .ensure_rrr_engine_registered()

  scenarios <- expand.grid(
    ar1 = c(FALSE, TRUE),
    nruns = c(1L, 2L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenarios))) {
    ar1 <- isTRUE(scenarios$ar1[[i]])
    nruns <- as.integer(scenarios$nruns[[i]])

    sim <- .simulate_rrr_parity_dataset(
      seed = 100L + i,
      nruns = nruns,
      ar1 = ar1,
      T_run = 60L,
      V = 40L
    )
    dset <- sim$dataset

    ar_opts <- if (ar1) list(struct = "ar1") else list(struct = "iid")
    fit_std <- fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      ar_options = ar_opts
    )
    ei <- fit_std$result$event_indices
    k_task <- length(ei)

    fit_rrr_full <- fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      ar_options = ar_opts,
      engine = "rrr_gls",
      engine_args = list(rank = k_task)
    )
    fit_rrr_red <- fmri_lm(
      onsets ~ hrf(condition),
      block = ~run,
      dataset = dset,
      ar_options = ar_opts,
      engine = "rrr_gls",
      engine_args = list(rank = 1L)
    )

    B_std <- .event_beta_matrix(fit_std)
    B_full <- .event_beta_matrix(fit_rrr_full)
    B_red <- .event_beta_matrix(fit_rrr_red)

    err_full <- mean(abs(B_std[ei, , drop = FALSE] - B_full[ei, , drop = FALSE]))
    err_red <- mean(abs(B_std[ei, , drop = FALSE] - B_red[ei, , drop = FALSE]))
    corr_full <- suppressWarnings(
      cor(
        as.numeric(B_std[ei, , drop = FALSE]),
        as.numeric(B_full[ei, , drop = FALSE])
      )
    )

    expect_true(is.finite(err_full))
    expect_true(is.finite(err_red))
    expect_true(is.finite(corr_full))

    # Single-run: near-exact parity. Multi-run: strong but not exact parity.
    if (nruns == 1L) {
      expect_lt(err_full, 1e-8)
      expect_gt(corr_full, 0.9999)
    } else {
      expect_lt(err_full, 0.08)
      expect_gt(corr_full, 0.95)
    }

    # Reduced-rank approximation should not outperform exact parity target.
    expect_gte(err_red + 1e-12, err_full)

    # Task-contrast estimate parity in full-rank setting.
    if (length(ei) >= 2L) {
      cvec <- c(1, -1)
      est_std <- as.numeric(cvec %*% B_std[ei[1:2], , drop = FALSE])
      est_full <- as.numeric(cvec %*% B_full[ei[1:2], , drop = FALSE])
      if (nruns == 1L) {
        expect_lt(max(abs(est_std - est_full)), 1e-8)
      } else {
        expect_gt(cor(est_std, est_full), 0.85)
        expect_lt(mean(abs(est_std - est_full)), 0.1)
      }
    }
  }
})

test_that("latent_sketch shows weak parity to exact task effects", {
  skip_on_cran()
  .ensure_rrr_engine_registered()

  sim <- .simulate_rrr_parity_dataset(
    seed = 301L,
    nruns = 1L,
    ar1 = FALSE,
    T_run = 80L,
    V = 64L,
    signal_sd = 0.6,
    noise_sd = 0.7
  )
  dset <- sim$dataset

  fit_std <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset
  )
  B_exact <- .event_beta_matrix(fit_std)
  ei <- fit_std$result$event_indices

  p <- ncol(B_exact)
  Tlen <- nrow(fmridataset::get_data_matrix(dset))
  low <- lowrank_control(
    time_sketch = list(method = "srht", m = min(8L * p, Tlen))
  )

  fit_sketch <- fmri_lm(
    onsets ~ hrf(condition),
    block = ~run,
    dataset = dset,
    engine = "latent_sketch",
    lowrank = low,
    ar_options = list(struct = "iid")
  )

  B_sketch <- fit_sketch$betas_fixed
  expect_equal(dim(B_sketch), dim(B_exact))
  expect_true(all(is.finite(B_sketch)))
  expect_true(all(is.finite(fit_sketch$sigma2)))

  corr <- suppressWarnings(
    cor(
      as.numeric(B_exact[ei, , drop = FALSE]),
      as.numeric(B_sketch[ei, , drop = FALSE])
    )
  )
  expect_true(is.finite(corr))
  expect_gt(corr, 0.05)
})
