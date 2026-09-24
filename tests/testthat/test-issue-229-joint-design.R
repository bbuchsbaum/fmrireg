make_issue_229_fixture <- function(equal_width = FALSE) {
  set.seed(229)
  n <- c(80L, 90L)
  sf <- fmrihrf::sampling_frame(n, TR = 2)
  ev <- data.frame(onset = rep(seq(8, 120, by = 16), 2),
                   run = rep(1:2, each = 8),
                   cond = factor(rep(c("a", "b"), 8)))
  con <- contrast_set(pair_contrast(~cond == "a", ~cond == "b", name = "a_b"))
  em <- event_model(onset ~ hrf(cond, contrasts = con), data = ev,
                    block = ~run, sampling_frame = sf)
  nu <- list(cbind(m = rnorm(n[1]), spike = c(1, rep(0, n[1] - 1))),
             cbind(m = rnorm(n[2])))
  if (equal_width) nu[[2]] <- cbind(nu[[2]], spike = c(0, 1, rep(0, n[2] - 2)))
  bm <- baseline_model(basis = "constant", sframe = sf, nuisance_list = nu)
  X <- cbind(as.matrix(design_matrix(em)), as.matrix(design_matrix(bm)))
  B <- matrix(rnorm(ncol(X) * 3), ncol(X), 3)
  Y <- X %*% B + matrix(rnorm(sum(n) * 3), sum(n), 3)
  ds <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  list(model = fmri_model(em, bm, ds), dataset = ds, X = X, Y = Y, n = n)
}

test_that("joint AR preserves global nuisance columns and matches direct GLS", {
  for (equal_width in c(FALSE, TRUE)) {
    f <- make_issue_229_fixture(equal_width)
    cfg <- fmri_lm_control(estimation = estimation_spec("joint"),
                           noise = noise_spec("ar1", exact_first = TRUE))
    phi <- 0.45
    pre <- prepare_chunkwise_matrices(f$model, f$dataset, cfg, phi_fixed = phi)
    # Independent stationary AR(1) precision, with no dependence across runs.
    Q <- as.matrix(Matrix::bdiag(lapply(f$n, function(n) {
      solve(toeplitz(phi^(0:(n - 1))) / (1 - phi^2))
    })))
    V <- solve(crossprod(f$X, Q %*% f$X))
    B <- V %*% crossprod(f$X, Q %*% f$Y)
    E <- f$Y - f$X %*% B
    df <- nrow(f$X) - qr(f$X)$rank
    ans <- process_chunk(f$Y, pre, cfg)
    expect_equal(ncol(pre$X_global), ncol(f$X))
    expect_equal(pre$proj_global$dfres, df)
    expect_equal(unname(ans$betas), unname(B), tolerance = 1e-9)
    expect_equal(as.numeric(ans$sigma2), colSums(E * (Q %*% E)) / df,
                 tolerance = 1e-9)
    expect_equal(as.numeric(pre$proj_global$XtXinv), as.numeric(V),
                 tolerance = 1e-9)
  }
})

test_that("public joint AR and robust fits retain global df and post-hoc covariance", {
  f <- make_issue_229_fixture()
  for (kind in c("ar1", "ar2", "robust", "ar_robust")) {
    cfg <- fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec(if (kind %in% c("robust")) "iid" else if (kind == "ar2") "ar2" else "ar1"),
      robust = robust_spec(if (kind %in% c("robust", "ar_robust")) "huber" else "none")
    )
    fit <- fmri_lm(f$model, dataset = f$dataset, control = cfg,
                   compute = compute_spec(voxel_chunks = 2))
    expect_equal(fit$result$rdf, nrow(f$X) - qr(f$X)$rank, info = kind)
    expect_equal(fit$result$cov.unscaled, fit$result$covariance_model_basis)
    post <- fit_contrasts(fit, list(a_b = c(1, -1)))$a_b
    inside <- fit$result$contrasts$data[[1]]
    expect_equal(as.numeric(post$estimate), as.numeric(inside$estimate), tolerance = 1e-9)
    expect_equal(as.numeric(post$se), as.numeric(inside$se), tolerance = 1e-9)
  }
})

test_that("file-backed joint AR handles unequal run nuisance widths", {
  f <- make_issue_229_fixture()
  td <- tempfile("issue-229-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  rows <- split(seq_len(sum(f$n)), rep(1:2, f$n))
  files <- vapply(1:2, function(r) {
    file <- file.path(td, paste0("run", r, ".nii"))
    neuroim2::write_vec(neuroim2::NeuroVec(
      array(t(f$Y[rows[[r]], ]), c(3, 1, 1, f$n[r])),
      neuroim2::NeuroSpace(c(3, 1, 1, f$n[r]), spacing = c(2, 2, 2))
    ), file)
    file
  }, character(1))
  mask <- file.path(td, "mask.nii")
  neuroim2::write_vol(neuroim2::NeuroVol(
    array(1, c(3, 1, 1)), neuroim2::NeuroSpace(c(3, 1, 1), spacing = c(2, 2, 2))
  ), mask)
  ds <- nifti_frame(files, mask = mask, TR = 2, run_length = f$n)
  cfg <- fmri_lm_control(estimation = estimation_spec("joint"), noise = noise_spec("ar1"))
  fit <- fmri_lm(f$model, dataset = ds, control = cfg,
                 compute = compute_spec(voxel_chunks = 2))
  memory <- fmri_lm(f$model, dataset = f$dataset, control = cfg,
                    compute = compute_spec(voxel_chunks = 2))
  expect_equal(fit$result$rdf, nrow(f$X) - qr(f$X)$rank)
  # NIfTI storage may round the double-precision synthetic input.
  expect_equal(as.numeric(coef(fit)), as.numeric(coef(memory)), tolerance = 1e-6)
  expect_equal(as.numeric(fit_contrasts(fit, list(a_b = c(1, -1)))$a_b$se),
               as.numeric(fit$result$contrasts$data[[1]]$se), tolerance = 1e-9)
})
