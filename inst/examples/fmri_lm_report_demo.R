#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fmrireg)
})

# If running from a source checkout before reinstalling fmrireg, load the
# local report implementation so the demo still works.
if (!("report" %in% getNamespaceExports("fmrireg")) && file.exists(file.path("R", "report.R"))) {
  source(file.path("R", "report.R"))
}

if (!exists("report", mode = "function")) {
  stop("report() is not available. Reinstall fmrireg or run this script from the package source tree.")
}

if (!requireNamespace("tinytable", quietly = TRUE)) {
  stop("tinytable is required to run this example script.")
}

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  quarto_bin <- "/Applications/quarto/bin/quarto"
}
if (!file.exists(quarto_bin)) {
  stop("Quarto CLI is required to run this example script.")
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else file.path(getwd(), "report-examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260222)

make_matrix_fit <- function() {
  sim <- simulate_simple_dataset(
    ncond = 2,
    nreps = 20,
    TR = 2,
    snr = 1.5,
    seed = 123
  )

  base <- sim$noisy[, -1, drop = FALSE]
  n_time <- nrow(base)
  n_vox <- 240L
  W <- matrix(rnorm(ncol(base) * n_vox), nrow = ncol(base), ncol = n_vox)
  Y <- base %*% W + matrix(rnorm(n_time * n_vox, sd = 0.35), nrow = n_time, ncol = n_vox)

  cond_vals <- factor(sim$conditions)
  if (nlevels(cond_vals) < 2L) {
    stop("Expected at least two conditions for contrast examples.")
  }
  levels(cond_vals)[seq_len(min(2L, nlevels(cond_vals)))] <- c("A", "B")[seq_len(min(2L, nlevels(cond_vals)))]

  etab <- data.frame(
    onset = sim$onsets,
    condition = cond_vals,
    run = 1L
  )

  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = 2,
    run_length = n_time,
    event_table = etab
  )

  con <- contrast_set(
    pair_contrast(~condition == "A", ~condition == "B", name = "A_vs_B")
  )

  fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~run,
    dataset = dset,
    durations = 0,
    progress = FALSE
  )
}

make_spatial_fit <- function() {
  dims <- c(18L, 18L, 10L)
  n_time <- 96L
  TR <- 2
  time_grid <- seq(0, by = TR, length.out = n_time)

  onsets <- seq(8, by = 6, length.out = 20)
  cond <- factor(rep(c("A", "B"), length.out = length(onsets)))
  etab <- data.frame(onset = onsets, condition = cond, run = 1L)

  regA <- fmrihrf::evaluate(
    fmrihrf::regressor(onsets = onsets[cond == "A"], hrf = fmrihrf::HRF_SPMG1, duration = 0, amplitude = 1),
    grid = time_grid
  )
  regB <- fmrihrf::evaluate(
    fmrihrf::regressor(onsets = onsets[cond == "B"], hrf = fmrihrf::HRF_SPMG1, duration = 0, amplitude = 1),
    grid = time_grid
  )

  gx <- array(rep(seq_len(dims[1]), each = dims[2] * dims[3]), dim = dims)
  gy <- array(rep(rep(seq_len(dims[2]), each = dims[3]), times = dims[1]), dim = dims)
  gz <- array(rep(seq_len(dims[3]), times = dims[1] * dims[2]), dim = dims)

  center1 <- c(6, 6, 5)
  center2 <- c(13, 13, 6)
  sphere1 <- (gx - center1[1])^2 + (gy - center1[2])^2 + (gz - center1[3])^2 <= 3^2
  sphere2 <- (gx - center2[1])^2 + (gy - center2[2])^2 + (gz - center2[3])^2 <= 3^2

  betaA <- array(0, dim = dims)
  betaB <- array(0, dim = dims)
  betaA[sphere1] <- 1.8
  betaB[sphere1] <- -0.4
  betaA[sphere2] <- -0.5
  betaB[sphere2] <- 1.9

  scans_arr <- array(0, dim = c(dims, n_time))
  for (tt in seq_len(n_time)) {
    signal_map <- betaA * regA[tt] + betaB * regB[tt]
    scans_arr[, , , tt] <- signal_map + array(rnorm(prod(dims), sd = 0.8), dim = dims)
  }

  scan_space <- neuroim2::NeuroSpace(dim = c(dims, n_time))
  mask_space <- neuroim2::NeuroSpace(dim = dims)
  scans <- list(neuroim2::NeuroVec(scans_arr, scan_space))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dim = dims), mask_space)

  dset <- fmridataset::fmri_mem_dataset(
    scans = scans,
    mask = mask,
    TR = TR,
    event_table = etab
  )

  con <- contrast_set(
    pair_contrast(~condition == "A", ~condition == "B", name = "A_vs_B"),
    pair_contrast(~condition == "B", ~condition == "A", name = "B_vs_A")
  )

  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = con),
    block = ~run,
    dataset = dset,
    durations = 0,
    strategy = "chunkwise",
    nchunks = 6L,
    progress = FALSE
  )

  bg_vol <- neuroim2::NeuroVol(apply(scans_arr, 1:3, mean), mask_space)

  synthetic_atlas <- list(
    atlas = {
      a <- array(0L, dim = dims)
      a[sphere1] <- 1L
      a[sphere2] <- 2L
      a
    },
    ids = c(1L, 2L),
    labels = c("Synthetic Region A", "Synthetic Region B")
  )

  list(fit = fit, bg = bg_vol, atlas = synthetic_atlas)
}

cat("Generating non-spatial report...\n")
fit_matrix <- make_matrix_fit()
matrix_pdf <- file.path(out_dir, "fmri_lm_report_matrix.pdf")
report(
  fit_matrix,
  output_file = matrix_pdf,
  title = "fMRI GLM Report (Matrix Dataset)",
  sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
  open = FALSE,
  quiet = TRUE
)
cat("  ->", matrix_pdf, "\n")

cat("Generating spatial report...\n")
sp <- make_spatial_fit()
spatial_pdf <- file.path(out_dir, "fmri_lm_report_spatial.pdf")
report(
  sp$fit,
  output_file = spatial_pdf,
  title = "fMRI GLM Report (Spatial Synthetic Dataset)",
  bg_vol = sp$bg,
  threshold = NULL,
  open = FALSE,
  quiet = TRUE
)
cat("  ->", spatial_pdf, "\n")

cat("Generating spatial + atlas-labeled report...\n")
atlas_pdf <- file.path(out_dir, "fmri_lm_report_spatial_atlas.pdf")
report(
  sp$fit,
  output_file = atlas_pdf,
  title = "fMRI GLM Report (Spatial + Synthetic Atlas Labels)",
  bg_vol = sp$bg,
  atlas = sp$atlas,
  threshold = NULL,
  open = FALSE,
  quiet = TRUE
)
cat("  ->", atlas_pdf, "\n")

cat("\nDone. Reports written to:", normalizePath(out_dir), "\n")
