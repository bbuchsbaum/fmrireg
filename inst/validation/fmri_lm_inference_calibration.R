# End-to-end calibration harness for fmri_lm inference.
#
# Every Monte Carlo draw regenerates the data and reruns AR estimation,
# whitening, residual covariance estimation, variance calculation, and
# studentization. This is intentionally not a conditional simulation from a
# covariance estimate held fixed across draws.

.fmri_lm_wilson_interval <- function(successes, trials, level = 0.99) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  phat <- successes / trials
  den <- 1 + z^2 / trials
  centre <- (phat + z^2 / (2 * trials)) / den
  half <- z * sqrt(phat * (1 - phat) / trials + z^2 / (4 * trials^2)) / den
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

run_fmri_lm_inference_calibration <- function(
    nsim = 500L,
    n = 120L,
    phi = 0.45,
    variance_method = c("model", "hac"),
    seed = 20260808L,
    absolute_tolerance = 0.025) {
  variance_method <- match.arg(variance_method)
  nsim <- as.integer(nsim)
  n <- as.integer(n)
  stopifnot(nsim >= 2L, n >= 60L, abs(phi) < 1)
  set.seed(seed)

  events <- data.frame(
    onset = seq(5, n - 25, length.out = 10),
    condition = factor(rep(c("A", "B"), 5)),
    run = 1L
  )
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    noise = if (variance_method == "model") noise_spec("ar1", iter_gls = 2L) else noise_spec("iid"),
    variance = if (variance_method == "model") {
      variance_spec("model", max_lag = 5L, df = "satterthwaite")
    } else {
      variance_spec("hac", max_lag = 5L, taper = "none",
                    df = "satterthwaite")
    }
  )
  compute <- compute_spec(voxel_chunks = 1L)

  draws <- matrix(NA_real_, nsim, 4L,
                  dimnames = list(NULL, c("estimate", "se", "p", "df")))
  for (i in seq_len(nsim)) {
    y <- as.numeric(stats::arima.sim(model = list(ar = phi), n = n))
    dset <- matrix_dataset(
      matrix(y, ncol = 1L), TR = 1, run_length = n,
      event_table = events
    )
    fit <- suppressWarnings(fmri_lm(
      onset ~ hrf(condition), block = ~run, dataset = dset,
      control = control, compute = compute
    ))
    payload <- fit$result$betas$data[[1L]]
    draws[i, ] <- c(
      payload$estimate[[1L]][1L, 1L],
      payload$se[[1L]][1L, 1L],
      payload$prob[[1L]][1L, 1L],
      fit$result$df$inference[1L]
    )
  }

  finite <- apply(draws, 1L, function(x) all(is.finite(x)))
  draws <- draws[finite, , drop = FALSE]
  if (nrow(draws) < 0.95 * nsim) {
    stop("More than 5% of calibration draws produced non-finite inference.",
         call. = FALSE)
  }
  q90 <- stats::qt(0.95, draws[, "df"])
  q95 <- stats::qt(0.975, draws[, "df"])
  metrics <- c(
    type1_05 = mean(draws[, "p"] < 0.05),
    coverage_90 = mean(abs(draws[, "estimate"]) <= q90 * draws[, "se"]),
    coverage_95 = mean(abs(draws[, "estimate"]) <= q95 * draws[, "se"])
  )
  targets <- c(type1_05 = 0.05, coverage_90 = 0.90, coverage_95 = 0.95)
  intervals <- t(vapply(seq_along(metrics), function(i) {
    .fmri_lm_wilson_interval(round(metrics[i] * nrow(draws)), nrow(draws))
  }, numeric(2)))
  variance_ratio <- mean(draws[, "se"]^2) / stats::var(draws[, "estimate"])

  table <- data.frame(
    metric = names(metrics), target = unname(targets), estimate = unname(metrics),
    lower_99 = intervals[, "lower"], upper_99 = intervals[, "upper"],
    pass = targets >= intervals[, "lower"] & targets <= intervals[, "upper"] &
      abs(metrics - targets) <= absolute_tolerance,
    row.names = NULL
  )

  structure(list(
    settings = list(nsim = nsim, n = n, phi = phi,
                    variance_method = variance_method, seed = seed),
    metrics = table,
    relative_variance = variance_ratio,
    relative_variance_pass = is.finite(variance_ratio) &&
      variance_ratio >= 0.85 && variance_ratio <= 1.15,
    draws = draws,
    pass = all(table$pass) && is.finite(variance_ratio) &&
      variance_ratio >= 0.85 && variance_ratio <= 1.15
  ), class = "fmri_lm_calibration")
}
