# Reproducible performance and retained-state guardrail for fmri_lm().
#
# Wall-clock timing is host-qualified, so the default slowdown threshold is
# deliberately broad. The hard contracts are numerical invariance across
# voxel partitions and the absence of transient inference residuals from the
# returned model-based fit.

run_fmri_lm_performance_guardrail <- function(
    n = 180L,
    nvox = 200L,
    chunks = c(1L, 4L),
    seed = 20260808L,
    max_fit_to_data_ratio = 12,
    max_chunk_slowdown = 6) {
  n <- as.integer(n)
  nvox <- as.integer(nvox)
  chunks <- as.integer(chunks)
  stopifnot(n >= 60L, nvox >= 10L, length(chunks) >= 2L,
            all(chunks >= 1L), max_fit_to_data_ratio > 0,
            max_chunk_slowdown > 0)
  set.seed(seed)

  events <- data.frame(
    onset = seq(5, n - 20, length.out = 12),
    condition = factor(rep(c("A", "B"), 6)),
    run = 1L
  )
  Y <- matrix(stats::rnorm(n * nvox), nrow = n, ncol = nvox)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    noise = noise_spec("ar1", iter_gls = 2L),
    robust = robust_spec("huber", max_iter = 2L),
    variance = variance_spec("model", df = "residual")
  )

  fits <- vector("list", length(chunks))
  elapsed <- numeric(length(chunks))
  for (i in seq_along(chunks)) {
    timing <- system.time({
      fits[[i]] <- suppressWarnings(fmri_lm(
        onset ~ hrf(condition), block = ~run, dataset = dset,
        control = control,
        compute = compute_spec(voxel_chunks = chunks[[i]])
      ))
    })
    elapsed[[i]] <- unname(timing[["elapsed"]])
  }

  reference <- coef(fits[[1L]])
  invariant <- vapply(fits[-1L], function(fit) {
    isTRUE(all.equal(coef(fit), reference, tolerance = 1e-10,
                     check.attributes = FALSE))
  }, logical(1))
  transient_state_absent <- vapply(fits, function(fit) {
    is.null(fit$result$inference_context) && is.null(fit$result$residuals)
  }, logical(1))
  size_ratio <- vapply(fits, function(fit) {
    as.numeric(utils::object.size(fit)) / as.numeric(utils::object.size(Y))
  }, numeric(1))
  baseline <- max(elapsed[[1L]], 1e-3)
  slowdown <- elapsed / baseline

  checks <- c(
    partition_invariant = all(invariant),
    transient_state_absent = all(transient_state_absent),
    retained_size_bounded = all(is.finite(size_ratio)) &&
      all(size_ratio <= max_fit_to_data_ratio),
    chunk_slowdown_bounded = all(is.finite(slowdown)) &&
      all(slowdown <= max_chunk_slowdown)
  )
  structure(list(
    settings = list(n = n, nvox = nvox, chunks = chunks, seed = seed),
    elapsed_seconds = stats::setNames(elapsed, paste0("chunks_", chunks)),
    slowdown = stats::setNames(slowdown, paste0("chunks_", chunks)),
    fit_to_data_size_ratio = stats::setNames(size_ratio, paste0("chunks_", chunks)),
    checks = checks,
    pass = all(checks)
  ), class = "fmri_lm_performance_guardrail")
}
