# LEGACY-GROUP-DELETE (next minor): whole file — tests the legacy readers read_h5_chunk() and the legacy fmri_meta.group_data_nifti() volume-reading path (read_nifti_data).
# Review 2026-10, stream B: group-data readers. HDF5 handle lifetime and
# failure reporting (B9); NIfTI volumes read once per subject (P-B1).

# --- B9 -------------------------------------------------------------------

h5_mock <- function(n_subj = 3L, n_vox = 6L, fail = integer(0), drop_se = integer(0)) {
  closed <- new.env()
  set.seed(99)
  vals <- lapply(seq_len(n_subj), function(i) {
    list(beta = stats::rnorm(n_vox), se = stats::runif(n_vox, 0.1, 0.5))
  })
  paths <- sprintf("sub-%02d.h5", seq_len(n_subj))
  handles <- lapply(seq_len(n_subj), function(i) {
    force(i)
    closed[[paths[i]]] <- 0L
    list(
      read_voxels = function(idx, labels) {
        if (i %in% fail) stop("corrupt dataset")
        out <- list(beta = vals[[i]]$beta[idx])
        if (!(i %in% drop_se)) out$se <- vals[[i]]$se[idx]
        out
      },
      close_all = function() closed[[paths[i]]] <- closed[[paths[i]]] + 1L
    )
  })
  names(handles) <- paths
  gd <- list(paths = paths, file_type = "by_stat", stat = c("beta", "se"), df = NULL)
  list(gd = gd, handles = handles, closed = closed, vals = vals, paths = paths)
}

test_that("B9: read_h5_chunk closes every subject's handle exactly once", {
  skip_if_not_installed("fmristore")
  mk <- h5_mock()
  testthat::local_mocked_bindings(
    read_labeled_vec = function(path) mk$handles[[path]],
    .package = "fmristore"
  )
  out <- fmrireg:::read_h5_chunk(mk$gd, voxel_indices = c(2L, 5L))
  closes <- vapply(mk$paths, function(p) mk$closed[[p]], integer(1))
  expect_equal(unname(closes), rep(1L, length(mk$paths)))
  expect_equal(out$beta, t(vapply(mk$vals, function(v) v$beta[c(2, 5)], numeric(2))))
  expect_equal(out$se, t(vapply(mk$vals, function(v) v$se[c(2, 5)], numeric(2))))
})

test_that("B9: a failed read is an error naming the file (and the handle is still closed)", {
  skip_if_not_installed("fmristore")
  mk <- h5_mock(fail = 2L)
  testthat::local_mocked_bindings(
    read_labeled_vec = function(path) mk$handles[[path]],
    .package = "fmristore"
  )
  expect_error(fmrireg:::read_h5_chunk(mk$gd, voxel_indices = 1:3), "sub-02.h5")
  expect_equal(mk$closed[["sub-02.h5"]], 1L)

  mk2 <- h5_mock(drop_se = 3L)
  testthat::local_mocked_bindings(
    read_labeled_vec = function(path) mk2$handles[[path]],
    .package = "fmristore"
  )
  expect_error(fmrireg:::read_h5_chunk(mk2$gd, voxel_indices = 1:3), "sub-03.h5")
})

test_that("B9: an unreadable file on disk is an error, not NA", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("hdf5r")
  bad <- tempfile(fileext = ".h5")
  writeLines("not an hdf5 file", bad)
  on.exit(unlink(bad), add = TRUE)
  gd <- list(paths = bad, file_type = "by_stat", stat = c("beta", "se"), df = NULL)
  expect_error(suppressWarnings(fmrireg:::read_h5_chunk(gd, voxel_indices = 1:2)),
               "Failed to read HDF5")
})

# --- P-B1 -----------------------------------------------------------------

mock_nifti_gd_io <- function(beta, se, covariates = NULL) {
  S <- nrow(beta)
  P <- ncol(beta)
  dims <- c(P, 1L, 1L)
  store <- new.env()
  bp <- paste0("beta_", seq_len(S))
  sp <- paste0("se_", seq_len(S))
  for (i in seq_len(S)) {
    store[[bp[i]]] <- array(beta[i, ], dims)
    store[[sp[i]]] <- array(se[i, ], dims)
  }
  gd <- structure(
    list(beta_paths = bp, se_paths = sp, var_paths = NULL, t_paths = NULL, df = NULL,
         subjects = paste0("s", seq_len(S)), covariates = covariates, mask = NULL,
         mask_data = NULL, dim = dims, format = "nifti", n_subjects = S, n_voxels = P),
    class = c("group_data_nifti", "group_data")
  )
  list(gd = gd, store = store)
}

test_that("P-B1: fmri_meta on NIfTI data reads each volume once, whatever the chunk size", {
  set.seed(5)
  S <- 7L
  P <- 23L
  beta <- matrix(stats::rnorm(S * P), S, P)
  se <- matrix(stats::runif(S * P, 0.2, 0.6), S, P)
  mk <- mock_nifti_gd_io(beta, se, covariates = data.frame(x = stats::rnorm(S)))
  reads <- new.env()
  reads$n <- 0L
  testthat::local_mocked_bindings(
    read_nifti_data = function(path, mask = NULL) {
      reads$n <- reads$n + 1L
      mk$store[[path]]
    }
  )
  fits <- lapply(c(5L, 23L, 1000L), function(cs) {
    reads$n <- 0L
    fit <- fmri_meta(mk$gd, formula = ~ x, method = "pm", chunk_size = cs,
                     verbose = FALSE, n_threads = 1)
    # one beta and one se volume per subject
    expect_equal(reads$n, 2L * S, info = paste("chunk_size", cs))
    fit
  })
  for (k in 2:3) {
    expect_equal(fits[[k]]$coefficients, fits[[1]]$coefficients)
    expect_equal(fits[[k]]$se, fits[[1]]$se)
    expect_equal(fits[[k]]$tau2, fits[[1]]$tau2)
    expect_equal(fits[[k]]$cov$tri, fits[[1]]$cov$tri)
  }
  # and the chunked result equals one direct fit of all voxels
  ref <- fmri_meta_fit(beta, se^2, stats::model.matrix(~ x, mk$gd$covariates),
                       method = "pm", n_threads = 1)
  expect_equal(unname(fits[[1]]$coefficients), unname(t(ref$beta)))
})
