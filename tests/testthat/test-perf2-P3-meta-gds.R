# Performance pass 2, stream P3-2 / P3-3: the GDS meta path computes the
# group dataset once, fits once (in feature chunks of chunk_size) and derives
# fit-time contrasts from that fit; t-combination is one vectorised call.

p3_gds_fixture <- function(S = 12L, P = 7L, seed = 11) {
  set.seed(seed)
  df <- expand.grid(subject = sprintf("s%02d", seq_len(S)), roi = sprintf("R%d", seq_len(P)),
                    contrast = "c1", stringsAsFactors = FALSE)
  sid <- as.integer(sub("s", "", df$subject))
  df$group <- factor(ifelse(sid <= S / 2, "A", "B"))
  df$x <- stats::rnorm(S)[sid]
  df$se <- stats::runif(nrow(df), 0.1, 0.4)
  df$beta <- 0.2 + 0.3 * (df$group == "B") +
    stats::rnorm(nrow(df), sd = sqrt(df$se^2 + 0.05))
  df$df <- 30 + (seq_len(nrow(df)) %% 5)
  df$t <- df$beta / df$se
  gd <- fmrireg:::.group_data_gds(
    df, format = "csv",
    effect_cols = c(beta = "beta", se = "se", t = "t", df = "df"),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
    covariate_cols = c("group", "x")
  )
  list(df = df, gd = gd, S = S, P = P)
}

p3_roi_data <- function(fx, roi) {
  d <- fx$df[fx$df$roi == roi, , drop = FALSE]
  d[order(d$subject), , drop = FALSE]
}

test_that("P3-2: GDS meta (ivw) matches metafor::rma for FE/DL/PM incl. contrasts and cov", {
  skip_if_not_installed("fmrigds")
  skip_if_not_installed("metafor")
  fx <- p3_gds_fixture()
  ctrl <- list(tol = 1e-13, threshold = 1e-13, maxiter = 10000)
  L <- rbind(a = c(0, 1, 0), b = c(0, 1, -1))
  for (m in c("fe", "dl", "pm")) {
    fit <- fmrireg:::.fmri_meta_group_data_gds_impl(
      fx$gd, formula = ~ group + x, method = m, weights = "ivw",
      contrasts = L, verbose = FALSE
    )
    rois <- rownames(fit$coefficients)
    expect_equal(sort(rois), sort(unique(fx$df$roi)))
    for (r in rois) {
      d <- p3_roi_data(fx, r)
      ref <- metafor::rma(yi = d$beta, vi = d$se^2, mods = ~ group + x, data = d,
                          method = toupper(m), control = ctrl)
      info <- paste(m, r)
      expect_equal(unname(fit$coefficients[r, ]), as.numeric(ref$beta), tolerance = 1e-6, info = info)
      expect_equal(unname(fit$se[r, ]), as.numeric(ref$se), tolerance = 1e-6, info = info)
      expect_equal(fit$tau2[match(r, rois)], ref$tau2, tolerance = 1e-6, info = info)
      expect_equal(fit$Q[match(r, rois)], ref$QE, tolerance = 1e-6, info = info)
      vc <- stats::vcov(ref)
      packed <- unlist(lapply(1:3, function(a) vc[a, a:3]))
      expect_equal(as.numeric(fit$cov$tri[, match(r, rois)]), as.numeric(packed),
                   tolerance = 1e-6, info = info)
      est <- unname(drop(L %*% ref$beta))
      se <- sqrt(diag(L %*% vc %*% t(L)))
      expect_equal(unname(fit$contrasts$estimate[match(r, rois), ]), est, tolerance = 1e-6, info = info)
      expect_equal(unname(fit$contrasts$se[match(r, rois), ]), unname(se), tolerance = 1e-6, info = info)
    }
  }
})

test_that("P3-2: GDS meta with equal / custom weights matches lm(weights =)", {
  skip_if_not_installed("fmrigds")
  fx <- p3_gds_fixture()
  wc <- seq(0.5, 2, length.out = fx$S)
  for (w in c("equal", "custom")) {
    fit <- fmrireg:::.fmri_meta_group_data_gds_impl(
      fx$gd, formula = ~ group, method = "fe", weights = w,
      weights_custom = if (w == "custom") wc else NULL,
      contrasts = c(groupB = 1), verbose = FALSE
    )
    for (r in rownames(fit$coefficients)) {
      d <- p3_roi_data(fx, r)
      ww <- if (w == "custom") wc else rep(1, fx$S)
      ref <- summary(stats::lm(beta ~ group, data = d, weights = ww))$coefficients
      i <- match(r, rownames(fit$coefficients))
      expect_equal(unname(fit$coefficients[i, ]), unname(ref[, 1]), tolerance = 1e-8)
      expect_equal(unname(fit$se[i, ]), unname(ref[, 2]), tolerance = 1e-8)
      expect_equal(unname(fit$contrasts$estimate[i, 1]), unname(ref[2, 1]), tolerance = 1e-8)
      expect_equal(unname(fit$contrasts$se[i, 1]), unname(ref[2, 2]), tolerance = 1e-8)
    }
  }
})

test_that("P3-2: results do not depend on chunk_size", {
  skip_if_not_installed("fmrigds")
  fx <- p3_gds_fixture()
  wm <- matrix(seq(0.5, 2, length.out = fx$S * fx$P), fx$S, fx$P)
  for (w in c("ivw", "custom")) {
    args <- list(fx$gd, formula = ~ group, method = "pm", weights = w,
                 weights_custom = if (w == "custom") wm else NULL,
                 contrasts = c(groupB = 1), verbose = FALSE)
    big <- do.call(fmrireg:::.fmri_meta_group_data_gds_impl, c(args, chunk_size = 1e4))
    small <- do.call(fmrireg:::.fmri_meta_group_data_gds_impl, c(args, chunk_size = 2))
    big$data <- small$data <- NULL
    expect_identical(small, big)
  }
})

test_that("P3-2 guard: one fmrigds::compute and one fit per call (incl. fit-time contrasts)", {
  skip_if_not_installed("fmrigds")
  fx <- p3_gds_fixture()
  n_compute <- 0L
  n_fit <- 0L
  real_compute <- fmrigds::compute
  real_fit <- fmrireg:::.meta_fit_core
  local_mocked_bindings(
    compute = function(...) { n_compute <<- n_compute + 1L; real_compute(...) },
    .package = "fmrigds"
  )
  local_mocked_bindings(
    .meta_fit_core = function(...) { n_fit <<- n_fit + 1L; real_fit(...) }
  )
  fit <- fmrireg:::.fmri_meta_group_data_gds_impl(
    fx$gd, formula = ~ group, method = "pm", weights = "ivw",
    contrasts = c(groupB = 1), chunk_size = 1e4, verbose = FALSE
  )
  expect_false(is.null(fit$contrasts))
  expect_equal(n_compute, 1L)
  expect_equal(n_fit, 1L)

  # chunk_size bounds the features per kernel call
  n_compute <- 0L; n_fit <- 0L
  fmrireg:::.fmri_meta_group_data_gds_impl(
    fx$gd, formula = ~ group, method = "pm", weights = "ivw",
    contrasts = c(groupB = 1), chunk_size = 3, verbose = FALSE
  )
  expect_equal(n_compute, 1L)
  expect_equal(n_fit, ceiling(fx$P / 3))
})

test_that("P3-3: vectorised t-combination equals per-feature combination and Stouffer closed form", {
  skip_if_not_installed("fmrigds")
  fx <- p3_gds_fixture()
  wc <- seq(0.5, 2, length.out = fx$S)
  wm <- matrix(seq(0.5, 2, length.out = fx$S * fx$P), fx$S, fx$P)
  cases <- list(
    list("stouffer", "equal", NULL), list("stouffer", "custom", wc),
    list("stouffer", "custom", wm), list("lancaster", "equal", NULL),
    list("lancaster", "custom", wm), list("fisher", "equal", NULL)
  )
  for (cs in cases) {
    fit <- fmrireg:::.fmri_meta_group_data_gds_impl(
      fx$gd, formula = ~ 1, combine = cs[[1]], weights = cs[[2]],
      weights_custom = cs[[3]], verbose = FALSE
    )
    rois <- rownames(fit$coefficients)
    for (r in rois) {
      d <- p3_roi_data(fx, r)
      j <- match(r, sort(unique(fx$df$roi)))
      wj <- if (is.null(cs[[3]])) NULL else if (is.matrix(cs[[3]])) cs[[3]][, j] else cs[[3]]
      per_feature <- fmrireg:::combine_t_statistics(
        matrix(d$t, ncol = 1), df = d$df, method = cs[[1]], weights = cs[[2]],
        weights_custom = wj
      )
      expect_equal(unname(fit$coefficients[r, 1]), per_feature, tolerance = 1e-12,
                   info = paste(cs[[1]], cs[[2]], r))
      if (cs[[1]] == "stouffer") {
        ww <- if (is.null(wj)) rep(1, fx$S) else wj
        zi <- stats::qnorm(stats::pt(d$t, d$df))
        expect_equal(unname(fit$coefficients[r, 1]), sum(ww * zi) / sqrt(sum(ww^2)),
                     tolerance = 1e-10)
      }
    }
  }
})

test_that("P3-3: combine_t_statistics accepts per-feature df and drops missing subjects", {
  set.seed(2)
  tm <- matrix(stats::rnorm(30), 6, 5)
  dfm <- matrix(10 + 1:30, 6, 5)
  tm[2, 3] <- NA
  dfm[4, 5] <- NA
  for (m in c("stouffer", "fisher", "lancaster")) {
    vec <- fmrireg:::combine_t_statistics(tm, df = dfm, method = m)
    per <- vapply(1:5, function(j) {
      fmrireg:::combine_t_statistics(matrix(tm[, j], ncol = 1), df = dfm[, j], method = m)
    }, numeric(1))
    expect_equal(vec, per, tolerance = 1e-14, info = m)
  }
  expect_error(fmrireg:::combine_t_statistics(tm, df = dfm[, 1:4]), "df")
})

test_that("P3-3 guard: one combine_t_statistics call for all features", {
  skip_if_not_installed("fmrigds")
  fx <- p3_gds_fixture()
  n_comb <- 0L
  real <- fmrireg:::combine_t_statistics
  local_mocked_bindings(
    combine_t_statistics = function(...) { n_comb <<- n_comb + 1L; real(...) }
  )
  fmrireg:::.fmri_meta_group_data_gds_impl(
    fx$gd, formula = ~ 1, combine = "stouffer", weights = "equal", verbose = FALSE
  )
  expect_equal(n_comb, 1L)
})
