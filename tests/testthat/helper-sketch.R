# Fixtures for the time-sketched (latent_sketch) engine.
#
# Sketch-and-solve needs more sketch rows than design columns, so its
# fixtures need more scans than the 8-scan .demo_matrix_dataset(), whose
# two-run design already has 10 columns (T < p).

# Two runs of `run_length` scans, two conditions, `nvox` voxels of iid N(0, 1)
# noise around `baseline`, plus condition effects drawn N(0, signal^2) per
# voxel when `signal` is non-zero.
sketch_matrix_dataset <- function(run_length = 30L, nvox = 2L, baseline = 100,
                                  signal = 0, seed = 1L) {
  withr::with_seed(seed, {
    TR <- 2
    nrun <- 2L
    onsets <- seq(4, run_length * TR - 12, by = 12)
    ev <- data.frame(
      onsets = rep(onsets, nrun),
      condition = factor(rep(rep(c("A", "B"), length.out = length(onsets)), nrun)),
      run = rep(seq_len(nrun), each = length(onsets))
    )
    Tlen <- nrun * run_length
    Y <- matrix(stats::rnorm(Tlen * nvox), Tlen, nvox) + baseline
    if (signal != 0) {
      sf <- fmridesign::sampling_frame(blocklens = rep(run_length, nrun), TR = TR)
      em <- fmridesign::event_model(onsets ~ hrf(condition), data = ev,
                                    block = ~run, sampling_frame = sf)
      Xe <- as.matrix(fmridesign::design_matrix(em))
      Y <- Y + Xe %*% matrix(signal * stats::rnorm(ncol(Xe) * nvox), ncol(Xe), nvox)
    }
    matrix_frame(Y, TR = TR, run_length = rep(run_length, nrun), event_table = ev)
  })
}
