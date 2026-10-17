# Group delegation to fmrigds, step P1 (data-raw/simplification-2026-10-plan.md,
# Phase 4 item 1). The gds path of fmri_meta() is pinned field by field
# against the unchanged compiled kernel (.meta_fit_core), which is in turn
# anchored to closed forms.

.p1_long <- function(S = 14L, P = 25L, K = 1L, seed = 11L) {
  set.seed(seed)
  subj <- sprintf("s%02d", seq_len(S))
  age <- round(runif(S, 20, 60))
  grp <- factor(rep(c("A", "B"), length.out = S))
  cons <- paste0("c", seq_len(K))
  df <- expand.grid(subject = subj, roi = sprintf("R%02d", seq_len(P)),
                    contrast = cons, stringsAsFactors = FALSE)
  i <- match(df$subject, subj)
  r <- match(df$roi, sprintf("R%02d", seq_len(P)))
  k <- match(df$contrast, cons)
  df$age <- age[i]
  df$group <- grp[i]
  df$se <- runif(nrow(df), 0.1, 0.4)
  tau <- rep(c(0, 0.1, 0.3, 0.6, 1), length.out = P)[r]
  df$beta <- 0.2 * r / P + 0.01 * (df$age - 40) + 0.3 * (df$group == "B") +
    0.5 * (k - 1) + rnorm(nrow(df), 0, sqrt(df$se^2 + tau^2))
  df$var <- df$se^2
  df
}
.p1_gds <- function(df, effect = c("se", "var")) {
  effect <- match.arg(effect)
  fmrireg::group_data(
    df, format = "csv",
    effect_cols = c(beta = "beta", stats::setNames(effect, effect)),
    subject_col = "subject", roi_col = "roi", contrast_col = "contrast",
    covariate_cols = c("age", "group")
  )
}
.p1_configs <- function() {
  list(
    fe_se    = list(effect = "se",  args = list(formula = ~ 1, method = "fe")),
    dl_se    = list(effect = "se",  args = list(formula = ~ 1, method = "dl")),
    pm_se    = list(effect = "se",  args = list(formula = ~ 1, method = "pm")),
    reml_se  = list(effect = "se",  args = list(formula = ~ 1, method = "reml")),
    pm_var   = list(effect = "var", args = list(formula = ~ 1, method = "pm")),
    dl_var_cov = list(effect = "var", args = list(formula = ~ 1 + age, method = "dl")),
    pm_cov2  = list(effect = "se",  args = list(formula = ~ 1 + age + group, method = "pm")),
    fe_cov   = list(effect = "se",  args = list(formula = ~ group, method = "fe")),
    huber_pm = list(effect = "se",  args = list(formula = ~ 1 + group, method = "pm", robust = "huber")),
    huber_fe = list(effect = "se",  args = list(formula = ~ 1, method = "fe", robust = "huber")),
    equal    = list(effect = "se",  args = list(formula = ~ 1 + group, weights = "equal")),
    equal_huber = list(effect = "se", args = list(formula = ~ 1, weights = "equal", robust = "huber")),
    custom_vec = list(effect = "se", args = list(formula = ~ 1 + age, weights = "custom",
                                                 weights_custom = seq(1, 2, length.out = 14))),
    custom_mat = list(effect = "var", args = list(formula = ~ 1, weights = "custom",
                                                  weights_custom = matrix(seq(0.5, 3, length.out = 14 * 25), 14, 25))),
    con_named = list(effect = "se", args = list(formula = ~ 1 + age + group, method = "pm",
                                                contrasts = c(groupB = 1, age = -2))),
    con_matrix = list(effect = "se", args = list(formula = ~ 1 + group, method = "dl",
                                                 contrasts = rbind(diff = c(0, 1), sum = c(1, 1)))),
    con_unnamed = list(effect = "var", args = list(formula = ~ 1 + age, method = "fe",
                                                   contrasts = c(1, 10))),
    con_equal = list(effect = "se", args = list(formula = ~ 1 + group, weights = "equal",
                                                contrasts = c("(Intercept)" = 1, groupB = 0.5)))
  )
}

# Reference: the direct kernel on the same subject x feature matrices.
.p1_reference <- function(df, effect, args, contrast = "c1") {
  d <- df[df$contrast == contrast, ]
  subj <- sort(unique(d$subject))
  rois <- sort(unique(d$roi))
  Y <- matrix(d$beta[match(paste(rep(subj, length(rois)), rep(rois, each = length(subj))),
                           paste(d$subject, d$roi))], length(subj), length(rois))
  V <- matrix(d$var[match(paste(rep(subj, length(rois)), rep(rois, each = length(subj))),
                          paste(d$subject, d$roi))], length(subj), length(rois))
  cov_df <- d[match(subj, d$subject), c("age", "group")]
  X <- stats::model.matrix(args$formula, cov_df)
  method <- suppressMessages(fmrireg:::.meta_resolve_method(args$method %||% "pm"))
  ref <- fmrireg:::.meta_fit_core(Y, V, X, method = method,
                                  robust = args$robust %||% "none",
                                  weights = args$weights %||% "ivw",
                                  weights_custom = args$weights_custom)
  terms <- colnames(X)
  I2 <- ref$I2_fe
  fill <- is.na(I2) & is.finite(ref$Q_fe) & is.finite(ref$df)
  I2[fill] <- pmax((ref$Q_fe[fill] - ref$df[fill]) / pmax(ref$Q_fe[fill], .Machine$double.eps), 0)
  out <- list(
    coefficients = matrix(t(ref$beta), dimnames = list(rois, terms), ncol = length(terms)),
    se = matrix(t(ref$se), dimnames = list(rois, terms), ncol = length(terms)),
    tau2 = as.numeric(ref$tau2), I2 = as.numeric(I2), Q = as.numeric(ref$Q_fe),
    Q_df = as.numeric(ref$df), X = X, method = method, terms = terms,
    cov_tri = ref$cov_tri, rois = rois, S = length(subj)
  )
  con <- args$contrasts
  if (!is.null(con)) {
    if (is.matrix(con)) {
      Cmat <- t(con); cn <- rownames(con)
    } else {
      w <- if (is.null(names(con))) con else {
        x <- stats::setNames(numeric(length(terms)), terms); x[names(con)] <- con; x
      }
      Cmat <- matrix(w, ncol = 1L); cn <- "contrast1"
    }
    cf <- fmrireg:::.meta_contrasts_from_tri(ref$beta, ref$cov_tri, Cmat)
    # a contrast matrix labels the columns; a single vector contrast does not
    dn <- if (is.matrix(con)) list(NULL, cn) else NULL
    out$contrasts <- list(names = cn,
                          estimate = matrix(t(cf$estimate), ncol = length(cn), dimnames = dn),
                          se = matrix(t(cf$se), ncol = length(cn), dimnames = dn),
                          z = matrix(t(cf$z), ncol = length(cn), dimnames = dn))
  }
  out
}

.p1_expect_fit <- function(fit, ref, args, info) {
  tol <- 1e-10
  expect_s3_class(fit, "fmri_meta")
  expect_setequal(names(fit), c("coefficients", "se", "tau2", "I2", "Q", "Q_df", "model",
                                "method", "robust", "weights", "data", "formula", "n_voxels",
                                "n_subjects", "roi_names", "cov",
                                if (!is.null(args$contrasts)) "contrasts"))
  expect_equal(fit[["coefficients"]], ref$coefficients, tolerance = tol, info = info)
  expect_equal(fit[["se"]], ref$se, tolerance = tol, info = info)
  expect_equal(fit[["tau2"]], ref$tau2, tolerance = tol, info = info)
  expect_equal(fit[["I2"]], ref$I2, tolerance = tol, info = info)
  expect_equal(fit[["Q"]], ref$Q, tolerance = tol, info = info)
  expect_equal(fit[["Q_df"]], ref$Q_df, tolerance = tol, info = info)
  expect_equal(unname(fit[["model"]][["X"]]), unname(ref$X), ignore_attr = TRUE, info = info)
  expect_identical(colnames(fit[["model"]][["X"]]), ref$terms, info = info)
  expect_identical(fit[["model"]][["formula"]], args$formula, info = info)
  expect_identical(fit[["formula"]], args$formula, info = info)
  expect_identical(fit[["method"]], ref$method, info = info)
  expect_identical(fit[["robust"]], args$robust %||% "none", info = info)
  expect_identical(fit[["weights"]], args$weights %||% "ivw", info = info)
  expect_equal(fit[["n_voxels"]], length(ref$rois), info = info)
  expect_equal(fit[["n_subjects"]], ref$S, info = info)
  expect_identical(fit[["roi_names"]], ref$rois, info = info)
  expect_identical(fit[["cov"]][["type"]], "tri", info = info)
  expect_identical(fit[["cov"]][["coef_names"]], ref$terms, info = info)
  expect_equal(fit[["cov"]][["tri"]], ref$cov_tri, tolerance = tol, info = info)
  if (!is.null(args$contrasts)) {
    expect_identical(fit[["contrasts"]][["names"]], ref$contrasts$names, info = info)
    for (f in c("estimate", "se", "z")) {
      expect_equal(fit[["contrasts"]][[f]], ref$contrasts[[f]], tolerance = tol,
                   info = paste(info, f))
    }
  }
}

test_that("gds fmri_meta() equals the direct meta kernel field by field (P1 pin)", {
  skip_if_not_installed("fmrigds")
  df <- .p1_long()
  for (nm in names(.p1_configs())) {
    cf <- .p1_configs()[[nm]]
    fit <- suppressMessages(do.call(
      fmri_meta, c(list(.p1_gds(df, cf$effect)), cf$args, list(verbose = FALSE))
    ))
    .p1_expect_fit(fit, .p1_reference(df, cf$effect, cf$args), cf$args, nm)
  }
})

test_that("the meta kernel matches closed forms (FE weighted mean, equal-weight lm)", {
  df <- .p1_long()
  d <- df[df$roi == "R07" & df$contrast == "c1", ]
  d <- d[order(d$subject), ]
  fe <- fmrireg:::.meta_fit_core(matrix(d$beta), matrix(d$var), matrix(1, nrow(d)), method = "fe")
  w <- 1 / d$var
  expect_equal(as.numeric(fe$beta), sum(w * d$beta) / sum(w), tolerance = 1e-12)
  expect_equal(as.numeric(fe$se), 1 / sqrt(sum(w)), tolerance = 1e-12)
  expect_equal(as.numeric(fe$Q_fe), sum(w * (d$beta - sum(w * d$beta) / sum(w))^2),
               tolerance = 1e-12)

  X <- stats::model.matrix(~ 1 + group, d)
  eq <- fmrireg:::.meta_fit_core(matrix(d$beta), matrix(d$var), X, weights = "equal")
  ls <- summary(stats::lm(beta ~ 1 + group, d))$coefficients
  expect_equal(as.numeric(eq$beta), unname(ls[, 1]), tolerance = 1e-10)
  expect_equal(as.numeric(eq$se), unname(ls[, 2]), tolerance = 1e-10)
})

test_that("reducer fmrireg:meta_reg returns exactly the kernel's numbers through fmrigds", {
  skip_if_not_installed("fmrigds")
  expect_false(is.null(fmrigds::get_reducer("fmrireg:meta_reg")))
  # registering again never replaces it
  before <- fmrigds::get_reducer("fmrireg:meta_reg")
  expect_false(fmrireg:::.register_fmrireg_meta_reg())
  expect_identical(fmrigds::get_reducer("fmrireg:meta_reg"), before)

  set.seed(1)
  S <- 12; P <- 7
  subj <- paste0("s", 1:S)
  b <- array(rnorm(P * S), c(P, S, 1), dimnames = list(paste0("v", 1:P), subj, "c1"))
  v <- array(runif(P * S, 0.05, 0.5), dim(b), dimnames(b))
  cd <- data.frame(age = rnorm(S), row.names = subj)
  g <- fmrigds::with_col_data(fmrigds::as_gds(list(beta = b, var = v)), cd)
  X <- stats::model.matrix(~ age, cd)
  for (cfg in list(list(method = "dl", robust = "huber", weights = "1/var"),
                   list(method = "pm", robust = "none", weights = "1/var"),
                   list(method = "fe", robust = "none", weights = "custom",
                        custom_weights = seq(1, 3, length.out = S)))) {
    r <- fmrigds::compute(fmrigds::reduce(
      g, method = "fmrireg:meta_reg", formula = ~ age, weights = cfg$weights,
      options = cfg[c("method", "robust", "custom_weights")[c(TRUE, TRUE, !is.null(cfg$custom_weights))]]
    ))
    ref <- fmrireg:::.meta_fit_core(t(b[, , 1]), t(v[, , 1]), X, method = cfg$method,
                                    robust = cfg$robust,
                                    weights = if (cfg$weights == "1/var") "ivw" else "custom",
                                    weights_custom = cfg$custom_weights)
    info <- paste(cfg$method, cfg$weights)
    for (k in 1:2) {
      tn <- colnames(X)[k]
      expect_identical(as.numeric(fmrigds::assay(r, paste0("coef:", tn))), as.numeric(ref$beta[k, ]), info = info)
      expect_identical(as.numeric(fmrigds::assay(r, paste0("se_coef:", tn))), as.numeric(ref$se[k, ]), info = info)
      expect_identical(as.numeric(fmrigds::assay(r, paste0("p_coef:", tn))),
                       2 * stats::pnorm(abs(as.numeric(ref$z[k, ])), lower.tail = FALSE), info = info)
    }
    expect_identical(as.numeric(fmrigds::assay(r, "tau2")), as.numeric(ref$tau2), info = info)
    expect_identical(as.numeric(fmrigds::assay(r, "Q")), as.numeric(ref$Q_fe), info = info)
    expect_identical(as.numeric(fmrigds::assay(r, "df_res")), as.numeric(ref$df), info = info)
    expect_identical(as.numeric(fmrigds::assay(r, "cov:(Intercept):age")), ref$cov_tri[2, ], info = info)
    expect_identical(fmrigds::coef_cov_tri(r, "fmrireg:meta_reg", contrast = 1)$cov_tri, ref$cov_tri,
                     info = info)
  }
})

test_that("gds fmri_meta() fits every subject-level contrast, each as if alone", {
  skip_if_not_installed("fmrigds")
  df <- .p1_long(K = 2L)
  args <- list(formula = ~ 1 + group, method = "pm", contrasts = c(groupB = 1))
  fits <- suppressMessages(do.call(fmri_meta, c(list(.p1_gds(df, "se")), args)))
  expect_s3_class(fits, "fmri_meta_list")
  expect_identical(names(fits), c("c1", "c2"))
  expect_output(print(fits), "2 subject-level contrasts")
  for (k in c("c1", "c2")) {
    .p1_expect_fit(fits[[k]], .p1_reference(df, "se", args, contrast = k), args, k)
    alone <- do.call(fmri_meta, c(list(.p1_gds(df[df$contrast == k, ], "se")), args))
    expect_equal(fits[[k]][["coefficients"]], alone[["coefficients"]], tolerance = 1e-12)
  }
  # the second contrast is a different analysis (its true intercept is 0.5 higher)
  expect_gt(mean(fits[["c2"]][["coefficients"]][, 1] - fits[["c1"]][["coefficients"]][, 1]), 0.3)
})
