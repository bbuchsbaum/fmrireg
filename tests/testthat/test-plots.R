skip_if_not_installed("ggplot2")

# layers are found by geom, not position, so adding a layer does not break tests
layer_by_geom <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
  ggplot2::layer_data(p, i[length(i)])
}

make_plot_fixture <- function(nvox = 12L, basis = "spmg1", extra_factor = FALSE,
                              partial_alias = FALSE, ar = 0) {
  withr::local_seed(3)
  TR <- 2
  nscan <- 80L
  sframe <- fmrihrf::sampling_frame(blocklens = c(nscan, nscan), TR = TR)
  events <- do.call(rbind, lapply(1:2, function(r) {
    ons <- seq(8, nscan * TR - 24, by = 10)
    data.frame(run = r, onset = ons,
               condition = factor(rep(c("faces", "houses"), length.out = length(ons))))
  }))
  events$rt_c <- as.numeric(scale(seq_len(nrow(events)) %% 5))
  # `task` duplicates `condition` exactly, which makes the design rank deficient
  events$task <- factor(ifelse(events$condition == "faces", "a", "b"))
  # `cue` spans the same events as `condition` with independent levels, so
  # exactly one of its columns is aliased
  events$cue <- factor(rep(c("left", "right", "right", "left"), length.out = nrow(events)))
  emod <- fmridesign::event_model(onset ~ hrf(condition, basis = basis) + hrf(rt_c),
                                  data = events, block = ~ run, sampling_frame = sframe)
  X <- as.matrix(fmridesign::design_matrix(emod))
  B <- matrix(0, ncol(X), nvox)
  B[1, ] <- 1.5
  E <- matrix(rnorm(nrow(X) * nvox), nrow(X))
  if (ar != 0) E <- apply(E, 2, function(e) as.numeric(stats::filter(e, ar, method = "recursive")))
  Y <- X %*% B + E + 100
  dset <- matrix_frame(Y, TR = TR, run_length = c(nscan, nscan), event_table = events)
  form <- if (extra_factor) {
    onset ~ hrf(condition, basis = basis) + hrf(rt_c) + hrf(task)
  } else if (partial_alias) {
    onset ~ hrf(condition, basis = basis) + hrf(cue)
  } else {
    onset ~ hrf(condition, basis = basis) + hrf(rt_c)
  }
  fit <- if (ar != 0) {
    suppressWarnings(fmri_lm(form, block = ~ run, dataset = dset,
                             control = fmri_lm_control(noise = noise_spec("ar1"))))
  } else {
    suppressWarnings(fmri_lm(form, block = ~ run, dataset = dset))
  }
  list(fit = fit, model = fit$model, Y = Y, events = events, sframe = sframe,
       nscan = nscan)
}

test_that("design_map puts the run boundary between runs and reports rank", {
  fx <- make_plot_fixture()
  p <- design_map(fx$model)
  expect_s3_class(p, "ggplot")
  hl <- layer_by_geom(p, "GeomHline")
  # scale_y_reverse() negates built positions
  expect_equal(abs(unique(hl$yintercept)), fx$nscan + 0.5)
  expect_match(p$labels$subtitle, "full rank")
  # run intercepts are drawn grey (NA fill) rather than as saturated blocks
  expect_true(anyNA(p$data$value[p$data$group == "intercept"]))
})

test_that(".column_info labels multi-basis columns by their true condition", {
  fx <- make_plot_fixture(basis = "spmg3")
  info <- .column_info(fx$model)
  ev <- info[info$kind == "event" & info$group == "condition", ]
  DM <- as.matrix(design_matrix(fx$model))
  tgrid <- fmrihrf::samples(fx$sframe, global = TRUE)
  truth <- do.call(cbind, lapply(c("faces", "houses"), function(cond) {
    e <- fx$events[fx$events$condition == cond, ]
    R <- as.matrix(fmrihrf::evaluate(
      fmrihrf::regressor(fmrihrf::global_onsets(fx$sframe, e$onset, e$run),
                         fmrihrf::HRF_SPMG3), tgrid))
    colnames(R) <- paste0(cond, " [b", 1:3, "]")
    R
  }))
  best <- apply(abs(stats::cor(DM[, ev$column], truth)), 1,
                function(r) colnames(truth)[which.max(r)])
  expect_equal(unname(ev$label), unname(best))
})

test_that("plot.fmri_model never joins lines across runs and drops intercepts", {
  fx <- make_plot_fixture()
  p <- plot(fx$model)
  expect_s3_class(p, "ggplot")
  d <- p$data
  boundary <- max(fmrihrf::samples(fx$sframe, global = TRUE)[seq_len(fx$nscan)])
  spans <- tapply(d$time, d$grp, function(t) any(t <= boundary) && any(t > boundary))
  expect_false(any(spans))
  expect_false("intercept" %in% levels(d$panel))
})

test_that("correlation_map blanks cross-run pairs and uses the run's own scans", {
  fx <- make_plot_fixture()
  p <- correlation_map(fx$model, half_matrix = FALSE)
  d <- p$data
  info <- .column_info(fx$model)
  lab_run <- stats::setNames(info$run, info$label)
  rr <- lab_run[as.character(d$row)]
  cr <- lab_run[as.character(d$col)]
  expect_false(any(!is.na(rr) & !is.na(cr) & rr != cr))
  expect_false(any(info$label[info$kind == "intercept"] %in% d$col))

  # an event column against a run-1 drift column: correlation over run 1 only
  DM <- as.matrix(design_matrix(fx$model))
  ev_col <- info$column[info$label == "faces"]
  dr <- info[info$kind == "drift" & info$run %in% 1L, ][1, ]
  rows <- seq_len(fx$nscan)
  expected <- stats::cor(DM[rows, ev_col], DM[rows, dr$column])
  got <- d$r[d$row == "faces" & d$col == dr$label]
  expect_equal(got, expected, tolerance = 1e-8)
})

test_that("correlation_map keeps constant columns visible and fails clearly when empty", {
  fx <- make_plot_fixture()
  DM <- as.matrix(design_matrix(fx$model))
  DM[, 2] <- 0
  info <- .column_info(fx$model, DM)
  p <- .correlation_map_common(DM, info = info, runs = c(fx$nscan, fx$nscan),
                               half_matrix = FALSE)
  expect_true(info$label[2] %in% as.character(p$data$col))
  expect_match(p$labels$caption, "constant column")

  expect_error(
    .correlation_map_common(matrix(1, 10, 3, dimnames = list(NULL, c("a", "b", "c")))),
    "constant"
  )
})

test_that("correlation_map works for a baseline_model alone", {
  fx <- make_plot_fixture()
  p <- correlation_map(fx$model$baseline_model)
  expect_s3_class(p, "ggplot")
  expect_gt(nrow(p$data), 0)
})

test_that("coefficient plot intervals are estimate +/- qt * se", {
  fx <- make_plot_fixture()
  p <- ggplot2::autoplot(fx$fit, voxel = 2, level = 0.9)
  d <- p$data
  q <- stats::qt(0.95, d$dfi)
  expect_equal(d$hi - d$est, q * d$se, tolerance = 1e-10)
  expect_equal(d$est - d$lo, q * d$se, tolerance = 1e-10)
})

test_that("t-distribution view builds and reports exceedance per tail", {
  fx <- make_plot_fixture()
  p <- ggplot2::autoplot(fx$fit)
  expect_silent(ggplot2::ggplot_build(p))
  tt <- as.matrix(stats(fx$fit))
  crit <- stats::qt(1 - 0.001 / 2, stats::median(fx$fit$result$df$inference))
  expected <- unname(.pct(colMeans(tt > crit)))
  shown <- unlist(lapply(p$layers, function(l) {
    if (inherits(l$geom, "GeomText") && !is.null(l$data$pos)) as.character(l$data$pos)
  }))
  expect_setequal(shown, expected)
})

test_that("a single aliased level is left out of the HRF view, not drawn as zero", {
  fx <- make_plot_fixture(partial_alias = TRUE)
  tt <- as.matrix(stats(fx$fit))
  aliased <- colnames(tt)[colSums(is.finite(tt)) == 0]
  expect_length(aliased, 1L)
  p <- ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1)
  expect_no_warning(ggplot2::ggplot_build(p))
  lab <- .column_info(fx$model)$condition[.column_info(fx$model)$column == aliased]
  expect_false(lab %in% as.character(p$data$condition))
  expect_match(p$labels$caption, "Not estimable")
  expect_no_warning(ggplot2::ggplot_build(ggplot2::autoplot(fx$fit, voxel = 1)))
  pl <- ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1, direct_labels = FALSE)
  expect_no_warning(ggplot2::ggplot_build(pl))
})

test_that("AR fits show prewhitened residuals that are close to white", {
  fx <- make_plot_fixture(nvox = 40L, ar = 0.5)
  p <- ggplot2::autoplot(fx$fit, type = "residuals")
  expect_setequal(levels(p$data$series), c("OLS", "prewhitened"))
  expect_match(p$labels$caption, "lag 1 by run")
  lag1 <- p$data$r[p$data$series == "prewhitened" & round(p$data$x) == 1]
  lag1_ols <- p$data$r[p$data$series == "OLS" & round(p$data$x) == 1]
  expect_gt(lag1_ols, 0.2)
  expect_lt(abs(lag1), 0.15)
})

test_that("single-basis HRF labels do not claim an estimated latency", {
  fx <- make_plot_fixture()
  p <- ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1)
  txt <- layer_by_geom(p, "GeomText")$label
  expect_false(any(grepl("[0-9] s$", txt)))
  fx3 <- make_plot_fixture(basis = "spmg3")
  p3 <- ggplot2::autoplot(fx3$fit, type = "hrf", voxel = 1)
  txt3 <- layer_by_geom(p3, "GeomText")$label
  expect_true(any(grepl("[0-9] s$", txt3)))
})

test_that("single-voxel HRF band matches lm() for a single-run OLS fit", {
  fx <- make_plot_fixture(basis = "spmg3")
  rows <- seq_len(fx$nscan)
  ev1 <- fx$events[fx$events$run == 1, ]
  d1 <- matrix_frame(fx$Y[rows, ], TR = 2, run_length = fx$nscan, event_table = ev1)
  f1 <- fmri_lm(onset ~ hrf(condition, basis = "spmg3"), block = ~ run, dataset = d1)
  sa <- seq(0, 20, by = 0.5)
  p <- ggplot2::autoplot(f1, type = "hrf", voxel = 1, sample_at = sa)
  X1 <- as.matrix(design_matrix(f1$model))
  lf <- stats::lm(fx$Y[rows, 1] ~ X1 - 1)
  V <- stats::vcov(lf)[1:3, 1:3]
  G <- as.matrix(fmrihrf::HRF_SPMG3(sa))
  half <- stats::qt(0.975, lf$df.residual) * sqrt(rowSums((G %*% V) * G))
  d <- p$data[p$data$condition == levels(p$data$condition)[1], ]
  expect_equal(drop(G %*% stats::coef(lf)[1:3]), d$mean, tolerance = 1e-6)
  expect_equal(half, (d$hi - d$lo) / 2, tolerance = 1e-6)
})

test_that("time course refits run by run and guards R2 for constant voxels", {
  fx <- make_plot_fixture()
  p <- ggplot2::autoplot(fx$fit, type = "timecourse", voxel = 3)
  d <- p$data
  expect_equal(d$value[d$series == "observed"], as.numeric(fx$Y[, 3]))
  rows <- seq_len(fx$nscan)
  X <- as.matrix(design_matrix(fx$model))
  keep <- colSums(X[rows, ] != 0) > 0
  ref <- stats::lm.fit(X[rows, keep], fx$Y[rows, 3])$fitted.values
  expect_equal(d$value[d$series == "fitted"][rows], unname(ref), tolerance = 1e-8)

  Yc <- fx$Y
  Yc[, 2] <- 100
  dc <- matrix_frame(Yc, TR = 2, run_length = c(fx$nscan, fx$nscan),
                     event_table = fx$events)
  fc <- fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = dc)
  pc <- ggplot2::autoplot(fc, type = "timecourse", voxel = 2)
  expect_match(pc$labels$subtitle, "n/a")
})

test_that("single-voxel residual autocorrelation matches stats::acf", {
  fx <- make_plot_fixture()
  expect_equal(dim(.acf_cols(matrix(stats::rnorm(50)), 5)), c(5L, 1L))
  p <- ggplot2::autoplot(fx$fit, type = "residuals", voxel = 3)
  r <- p$data$r
  expect_equal(length(r), 10L)
  expect_gt(length(unique(round(r, 6))), 1L)
  # reference: per-run OLS residuals, autocorrelation pooled by run length
  X <- as.matrix(design_matrix(fx$model))
  ref <- rowMeans(sapply(1:2, function(run) {
    rows <- (run - 1) * fx$nscan + seq_len(fx$nscan)
    keep <- colSums(X[rows, ] != 0) > 0
    e <- stats::lm.fit(X[rows, keep], fx$Y[rows, 3])$residuals
    stats::acf(e, lag.max = 3, plot = FALSE)$acf[2:4]
  }))
  expect_equal(r[1:3], ref, tolerance = 1e-8)
})

test_that("every view builds silently for a rank-deficient fit", {
  fx <- make_plot_fixture(extra_factor = TRUE)
  expect_match(design_map(fx$model)$labels$subtitle, "RANK DEFICIENT")
  cap <- correlation_map(fx$model)$labels$caption
  expect_match(cap, "Inf \\(aliased\\)")
  expect_false(grepl("-[0-9]", cap))
  for (p in list(ggplot2::autoplot(fx$fit),
                 ggplot2::autoplot(fx$fit, voxel = 1),
                 ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1),
                 ggplot2::autoplot(fx$fit, type = "timecourse", voxel = 1),
                 ggplot2::autoplot(fx$fit, type = "residuals", voxel = 1:4))) {
    expect_no_warning(expect_no_error(ggplot2::ggplot_build(p)))
  }
})

test_that("fit plots validate voxel and level", {
  fx <- make_plot_fixture()
  expect_error(ggplot2::autoplot(fx$fit, voxel = NA_real_), "integer indices")
  expect_error(ggplot2::autoplot(fx$fit, voxel = 2.5), "integer indices")
  expect_error(ggplot2::autoplot(fx$fit, level = 2), "between 0 and 1")
  expect_error(ggplot2::autoplot(fx$fit, type = "hrf", direct_labels = NA), "TRUE or FALSE")
  expect_error(ggplot2::autoplot(fx$fit, type = "timecourse", voxel = 1:2), "single index")
})

test_that("plots draw on the pdf device without glyph warnings", {
  fx <- make_plot_fixture()
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  fx3 <- make_plot_fixture(basis = "spmg3")
  reg <- fmrihrf::regressor(onsets = c(10, 40, 45), hrf = fmrihrf::HRF_SPMG3)
  expect_no_warning({
    print(design_map(fx$model))
    print(correlation_map(fx$model))
    print(plot(fx$model))
    print(ggplot2::autoplot(reg))
    print(ggplot2::autoplot(fx$fit))
    print(ggplot2::autoplot(fx$fit, voxel = 1))
    print(ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1))
    print(ggplot2::autoplot(fx3$fit, type = "hrf", voxel = 1:4))
    print(ggplot2::autoplot(fx$fit, type = "timecourse", voxel = 1))
    print(ggplot2::autoplot(fx$fit, type = "residuals"))
  })
})

test_that("options(fmrireg.plot_theme) replaces the base theme", {
  fx <- make_plot_fixture()
  withr::local_options(fmrireg.plot_theme = ggplot2::theme_bw())
  p <- design_map(fx$model)
  expect_equal(p$theme$panel.border$colour, ggplot2::theme_bw()$panel.border$colour)
  pf <- ggplot2::autoplot(fx$fit, type = "hrf", voxel = 1)
  expect_equal(pf$theme$panel.border$colour, ggplot2::theme_bw()$panel.border$colour)
})

test_that("residual autocorrelation band is centred on the design-induced bias", {
  fx <- make_plot_fixture()
  withr::local_seed(21)
  n <- nrow(fx$Y)
  Yw <- matrix(stats::rnorm(n * 150), n)
  X <- as.matrix(design_matrix(fx$model))
  run_id <- rep(1:2, each = fx$nscan)
  sr <- .residual_acf(Yw, X, run_id, runwise = TRUE, lag_max = 5L)
  # least-squares residuals of white noise are negatively autocorrelated
  expect_lt(sr$expect[1], 0)
  expect_equal(rowMeans(sr$acf), sr$expect, tolerance = 0.02)
  # and the 95% band flags close to 5% of white-noise voxels at lag 1
  half <- stats::qnorm(0.975) * sqrt(diag(sr$cov))
  flag <- mean(abs(sr$acf[1, ] - sr$expect[1]) > half[1])
  expect_lt(flag, 0.12)
  # the exact null variance matches the spread of white-noise autocorrelations
  ratio <- apply(sr$acf - sr$expect, 1, stats::var) / diag(sr$cov)
  expect_true(all(ratio > 0.6 & ratio < 1.5))
})

test_that("every fmri_lm view uses the current accessor API without warnings", {
  # Deprecation warnings must be visible for this test to mean anything.
  withr::local_options(fmrireg.suppress_deprecation = FALSE)
  # 5 voxels x 2 event regressors: non-square, so a transposed voxels x terms
  # matrix would change dimensions or index the wrong cells.
  withr::local_seed(21)
  ev <- data.frame(onset = seq(10, 180, by = 12), run = 1,
                   condition = factor(rep(c("a", "b"), length.out = 15)))
  dset <- matrix_frame(matrix(rnorm(100 * 5), 100, 5), TR = 2, run_length = 100,
                       event_table = ev)
  con <- pair_contrast(~ condition == "a", ~ condition == "b", name = "a_vs_b")
  fit <- fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~ run, dataset = dset)
  expect_identical(dim(coef(fit)), c(5L, 2L))

  views <- list(
    list(), list(voxel = 2), list(type = "betas", voxel = 1:3),
    list(type = "contrasts"), list(type = "contrasts", voxel = 4),
    list(type = "hrf"), list(type = "hrf", voxel = 5),
    list(type = "timecourse", voxel = 3), list(type = "residuals"),
    list(type = "residuals", voxel = 1:2)
  )
  for (args in views) {
    expect_no_warning(p <- do.call(ggplot2::autoplot, c(list(fit), args)))
    expect_no_warning(ggplot2::ggplot_build(p))
  }
  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  expect_no_warning(plot(fit, type = "betas"))
  grDevices::dev.off()

  # The coefficient plot shows coef() for the chosen voxel, in coef()'s order.
  d <- ggplot2::autoplot(fit, voxel = 4)$data
  expect_equal(d$est, unname(coef(fit)[4, ]), tolerance = 1e-12)
  expect_equal(d$se, unname(as.matrix(standard_error(fit))[4, ]), tolerance = 1e-12)
  dc <- ggplot2::autoplot(fit, type = "contrasts", voxel = 4)$data
  expect_equal(dc$est, unname(as.matrix(coef(fit, type = "contrasts"))[4, "a_vs_b"]),
               tolerance = 1e-12)

  # "estimates" is kept as a silent synonym of "betas" for the plot type.
  expect_no_warning(old <- ggplot2::autoplot(fit, type = "estimates", voxel = 4))
  expect_equal(old$data, d)
})

test_that("coefficient plots index the non-square demo fit by voxel", {
  withr::local_options(fmrireg.suppress_deprecation = FALSE)
  fit <- suppressWarnings(.demo_fmri_lm())
  expect_identical(dim(coef(fit)), c(3L, 2L))
  for (v in 1:3) {
    expect_no_warning(d <- ggplot2::autoplot(fit, voxel = v)$data)
    expect_equal(d$est, unname(coef(fit)[v, ]), tolerance = 1e-12)
  }
  expect_no_warning(ggplot2::ggplot_build(ggplot2::autoplot(fit)))
})
