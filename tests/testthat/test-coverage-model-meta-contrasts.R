# Coverage for fmri_model methods, contrast API, metafuns, and df helpers.

tiny_fmri_model <- function(n_time = 80L, n_vox = 6L, n_events = 8L) {
  set.seed(11)
  run_len <- n_time / 2
  onsets <- c(
    seq(6, run_len - 8, length.out = n_events / 2),
    seq(6, run_len - 8, length.out = n_events / 2)
  )
  etab <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = n_events)),
    run = rep(1:2, each = n_events / 2)
  )
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  dset <- matrix_dataset(Y, TR = 1, run_length = c(n_time / 2, n_time / 2), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition),
    data = etab,
    block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 2, sframe = dset$sampling_frame)
  list(model = fmri_model(emod, bmod, dset), dataset = dset, events = etab)
}

test_that("fmri_model accessors, print, plot, design_map, correlation_map", {
  fx <- tiny_fmri_model()
  fmod <- fx$model

  expect_s3_class(fmod, "fmri_model")
  dm <- design_matrix(fmod)
  expect_equal(nrow(dm), 80L)
  expect_true(ncol(dm) >= 2)

  expect_true(length(terms(fmod)) >= 2)
  expect_true(length(event_terms(fmod)) >= 1)
  expect_true(length(baseline_terms(fmod)) >= 1)
  expect_equal(sum(blocklens(fmod)), 80L)
  expect_true(length(conditions(fmod)) >= 2)

  cel <- cells(fmod)
  expect_s3_class(cel, "tbl_df")

  out <- capture.output(print(fmod))
  expect_true(any(grepl("fMRI Model", out)))

  skip_if_not_installed("cowplot")
  p <- plot(fmod)
  expect_true(!is.null(p))

  cm <- correlation_map(fmod, half_matrix = TRUE)
  expect_s3_class(cm, "ggplot")

  dmap <- design_map(fmod, fill_midpoint = 0, rotate_x_text = FALSE)
  expect_s3_class(dmap, "ggplot")
  dmap2 <- design_map(fmod, fill_midpoint = NULL, block_separators = FALSE)
  expect_s3_class(dmap2, "ggplot")

  expect_error(fmrireg:::prediction_matrix(fmod), "not implemented")
  expect_error(fmrireg:::design_env.fmri_model(fmod), "Not implemented")

  # create_fmri_model path with explicit and NULL baseline
  created <- create_fmri_model(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = fx$dataset,
    durations = 0
  )
  expect_s3_class(created, "fmri_model")

  created2 <- create_fmri_model(
    onset ~ hrf(condition),
    block = ~ run,
    baseline_model = baseline_model(basis = "constant", sframe = fx$dataset$sampling_frame),
    dataset = fx$dataset
  )
  expect_s3_class(created2, "fmri_model")

  expect_error(
    create_fmri_model("not-a-formula", block = ~ run, dataset = fx$dataset),
    "formula"
  )
})

test_that("compute_lm_contrasts covers name-based t/F and validation errors", {
  set.seed(2)
  p <- 4L
  V <- 7L
  B <- matrix(rnorm(p * V), p, V)
  rownames(B) <- paste0("b", 1:p)
  XtXinv <- diag(p)
  colnames(XtXinv) <- rownames(XtXinv) <- rownames(B)
  sigma <- rep(1, V)

  stacked <- compute_lm_contrasts(
    B, XtXinv, df = 50,
    sigma = sigma,
    t_contrasts = list(A_vs_B = c(b1 = 1, b2 = -1)),
    f_contrasts = list(main = diag(2)),
    columns = rownames(B),
    output = "stacked"
  )
  expect_s3_class(stacked, "tbl_df")
  expect_true(nrow(stacked) >= 1)

  listed <- compute_lm_contrasts(
    B, XtXinv, df = 50,
    sigma2 = sigma^2,
    contrasts = list(t1 = c(1, -1, 0, 0), F1 = diag(2)),
    output = "list"
  )
  expect_true(is.list(listed))

  expect_error(compute_lm_contrasts(1:3, XtXinv, 10, sigma = 1), "matrix")
  expect_error(compute_lm_contrasts(B, XtXinv, 10), "sigma")
  expect_error(
    compute_lm_contrasts(B, XtXinv, 10, sigma = 1:2),
    "Length of sigma"
  )
})

test_that("metafuns stouffer/fixef/lancaster and meta_contrasts wiring", {
  set.seed(9)
  pval <- matrix(c(0.01, 0.2, 0.05, 0.4, 0, 0.1), nrow = 2)
  se <- matrix(c(0.2, 0.3, 0.25, 0.4, 0.2, 0.35), nrow = 2)
  st <- fmrireg:::meta_stouffer(pval, se)
  expect_equal(names(st), c("estimate", "se", "stat", "prob", "stat_type"))
  expect_equal(length(st$estimate), 2)

  se_m <- matrix(c(0.2, 0.25, 0.3, 0.35), 2, 2)
  beta_m <- matrix(c(0.5, 0.4, 0.6, 0.2), 2, 2)
  fe_inv <- fmrireg:::do_fixef(se_m, beta_m, "inv_var")
  expect_s3_class(fe_inv, "tbl_df")
  fe_eq <- fmrireg:::do_fixef(se_m, beta_m, "equal")
  expect_equal(nrow(fe_eq), 2)

  # Build contrast table for meta_fixef / meta_Fcontrasts
  mk_row <- function(name, type, estimate, se) {
    tibble::tibble(
      type = type,
      name = name,
      conmat = list(matrix(c(1, -1), 1)),
      colind = list(1:2),
      data = list(tibble::tibble(estimate = estimate, se = se, prob = c(0.05, 0.1)))
    )
  }
  ctab <- dplyr::bind_rows(
    mk_row("c1", "contrast", c(0.5, 0.4), c(0.2, 0.25)),
    mk_row("c1", "contrast", c(0.3, 0.2), c(0.22, 0.3))
  )
  out_fe <- fmrireg:::meta_fixef(ctab, weighting = "inv_var")
  expect_equal(out_fe$name[1], "c1")

  ftab <- dplyr::bind_rows(
    mk_row("f1", "Fcontrast", c(0.5, 0.4), c(1, 1)),
    mk_row("f1", "Fcontrast", c(0.3, 0.2), c(1, 1))
  )
  # meta_Fcontrasts expects data$prob
  ftab$data <- lapply(ftab$data, function(d) {
    d$prob <- c(0.04, 0.2)
    d
  })
  out_f <- fmrireg:::meta_Fcontrasts(ftab)
  expect_equal(out_f$stat_type[1], "meta_zfstat")

  # combine_t_statistics lancaster + warning paths
  tmat <- matrix(c(1.5, 2.0, 0.5, -1.2, 1.1, 0.8), nrow = 3)
  z_l <- fmrireg:::combine_t_statistics(
    tmat, df = 20, method = "lancaster", weights = "equal"
  )
  expect_equal(length(z_l), ncol(tmat))

  expect_warning(
    fmrireg:::combine_t_statistics(
      tmat, df = 20, method = "stouffer", weights = "ivw"
    ),
    "equal weights"
  )
  expect_warning(
    fmrireg:::combine_t_statistics(
      tmat, df = 20, method = "stouffer", weights = "custom"
    ),
    "equal weights"
  )
  expect_warning(
    fmrireg:::combine_t_statistics(
      tmat, df = 20, method = "fisher", weights = "ivw"
    ),
    "Weighted Fisher"
  )

  z_ivw <- fmrireg:::combine_t_statistics(
    tmat, df = 20, method = "stouffer", weights = "ivw",
    se_mat = matrix(0.3, 3, 2)
  )
  expect_equal(length(z_ivw), 2)

  expect_error(
    fmrireg:::combine_t_statistics(tmat, df = c(1, 2), method = "stouffer"),
    "Length of 'df'"
  )

  # Empty meta_contrasts
  empty <- fmrireg:::.meta_contrasts_impl(list(), weighting = "inv_var")
  expect_equal(nrow(empty), 1L)
})

test_that("effective df helpers cover AR/robust branches", {
  X <- cbind(1, rnorm(40), rnorm(40))
  fit_ols <- list()
  expect_equal(fmrireg:::compute_effective_df(fit_ols, X, "ols"), 40 - 3)

  fit_ar <- list(ar_coef = 0.4)
  df_ar <- fmrireg:::compute_effective_df(fit_ar, X, "ar")
  expect_true(df_ar < 40 - 3 || is.finite(df_ar))

  expect_equal(fmrireg:::compute_effective_df_ar(40, 3, NULL), 37)

  w <- runif(40, 0.4, 1)
  df_r <- fmrireg:::compute_effective_df_robust(X, w, solve(crossprod(X)))
  expect_true(df_r < 37)
  expect_equal(fmrireg:::compute_effective_df_robust(X, rep(1, 40), NULL), 37)

  fit_ar_r <- list(ar_coef = 0.3, robust_weights = w, XtXinv = solve(crossprod(X)))
  df_ar_r <- fmrireg:::compute_effective_df(fit_ar_r, X, "ar_robust")
  expect_true(is.finite(df_ar_r))
})

test_that("paired_diff_block and flip_sign cover ttest helpers", {
  set.seed(4)
  S <- 6
  P <- 5
  blkA <- list(
    Y = matrix(rnorm(S * P), S, P),
    V = matrix(runif(S * P, 0.05, 0.2), S, P),
    meta = list(subjects = paste0("s", 1:S), contrast = "A"),
    index = 1:P,
    covars = NULL,
    feature = "voxel"
  )
  blkB <- blkA
  blkB$Y <- matrix(rnorm(S * P), S, P)
  blkB$V <- matrix(runif(S * P, 0.05, 0.2), S, P)
  blkB$meta$contrast <- "B"

  diff <- paired_diff_block(blkA, blkB, rho = 0.2)
  expect_equal(dim(diff$Y), c(S, P))
  expect_equal(dim(diff$V), c(S, P))

  diff2 <- paired_diff_block(blkA, blkB, rho = rep(0.1, S))
  expect_equal(dim(diff2$V), c(S, P))

  expect_error(paired_diff_block(list(Y = 1), blkB), "Y matrices")
  expect_error(
    paired_diff_block(blkA, within(blkB, Y <- Y[, 1:2, drop = FALSE])),
    "identical dimensions"
  )

  # flip_sign on a synthetic fit object (flips by rownames)
  fit <- list(
    beta = matrix(c(1, -2, 3), 3, 1, dimnames = list(c("a", "b", "c"), NULL)),
    t = matrix(c(1.5, -2.5, 0.5), 3, 1, dimnames = list(c("a", "b", "c"), NULL)),
    z = matrix(c(1.2, -2.1, 0.4), 3, 1, dimnames = list(c("a", "b", "c"), NULL))
  )
  flipped <- flip_sign(fit, coef = "b")
  expect_equal(unname(flipped$beta["b", 1]), 2)
  expect_equal(unname(flipped$t["b", 1]), 2.5)
  flipped_all <- flip_sign(fit)
  expect_equal(unname(flipped_all$beta["a", 1]), -1)

  # weight validation helpers
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(NULL, 3, 4),
    "weights_custom"
  )
  w <- fmrireg:::.fmri_ttest_validate_weights_custom(c(1, 2, 3), 3, 4)
  expect_equal(w, c(1, 2, 3))
  expect_error(
    fmrireg:::.fmri_ttest_validate_weights_custom(c(1, -1, 2), 3, 4),
    "positive"
  )
  W <- fmrireg:::.fmri_ttest_validate_weights_custom(matrix(1, 3, 4), 3, 4)
  expect_equal(dim(W), c(3, 4))
})

test_that("group_data_from_fmrilm validates inputs", {
  etab <- data.frame(onset = c(5, 25), condition = factor(c("A", "B")), run = 1L)
  # matrix_dataset lacks NeuroSpace mask — expect spatial reconstruction error
  Y <- matrix(rnorm(50 * 4), 50, 4)
  dset <- matrix_dataset(Y, TR = 1, run_length = 50, event_table = etab)
  fit <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dset)

  expect_error(group_data_from_fmrilm("not-a-list"), "must be a list")
  expect_error(group_data_from_fmrilm(list(1)), "fmri_lm")

  # Spatial metadata failure is expected for matrix_dataset
  expect_error(
    suppressWarnings(group_data_from_fmrilm(list(fit, fit))),
    "Spatial metadata|fmri_dataset|neuroim2|mask"
  )
})

test_that("spatial mask helpers reject non-dataset inputs", {
  expect_error(
    fmrireg:::.fmri_dataset_mask_space(list(), "test"),
    "fmri_dataset"
  )
  expect_null(fmrireg:::.fmri_try_space(NULL))
  expect_null(fmrireg:::.fmri_try_space("path.nii"))
})
