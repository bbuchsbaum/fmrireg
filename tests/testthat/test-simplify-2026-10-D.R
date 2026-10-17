# Regression tests for the 2026-10 simplification, stream D (group analysis).

.d_demo_gds <- function() {
  gd_csv <- suppressWarnings(fmrireg:::.demo_group_data_csv())
  fmrireg::group_data(
    gd_csv$data, format = "csv",
    effect_cols = c(beta = "beta", se = "se"),
    subject_col = "subject",
    roi_col = "roi",
    contrast_col = "contrast",
    covariate_cols = "age"
  )
}

.d_registry_snapshot <- function() {
  red <- c("meta:fe", "meta:re", "meta:fe_reg", "meta:re_reg",
           "meta:re:dl", "meta:re:pm", "meta:re:reml",
           "frg:meta:fe", "frg:meta:re:dl", "frg:meta:re:pm", "frg:meta:re:reml")
  ph <- c("fdr:bh", "fdr:by", "fdr:spatial")
  list(
    reducers = fmrigds::list_reducers(),
    posthoc = fmrigds::list_posthoc(),
    reducer_fun = lapply(stats::setNames(red, red), function(n) {
      r <- fmrigds::get_reducer(n)
      if (is.null(r)) NULL else list(provides = r$provides, fun = r$fun)
    }),
    posthoc_fun = lapply(stats::setNames(ph, ph), function(n) {
      fmrigds::get_posthoc(n)$fun
    })
  )
}

test_that("gds fmri_meta()/fmri_ttest() leave fmrigds's reducer registry alone (B1)", {
  skip_if_not_installed("fmrigds")
  gd <- .d_demo_gds()
  before <- .d_registry_snapshot()
  # fmrigds's own meta:fe provides the *_g outputs
  expect_setequal(before$reducer_fun[["meta:fe"]]$provides,
                  c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"))
  fit <- fmri_meta(gd, formula = ~ 1, method = "pm", verbose = FALSE)
  tt <- fmri_ttest(gd, formula = ~ 1, engine = "meta")
  expect_s3_class(fit, "fmri_meta")
  after <- .d_registry_snapshot()
  expect_identical(after$reducers, before$reducers)
  expect_identical(after$posthoc, before$posthoc)
  for (nm in names(before$reducer_fun)) {
    expect_identical(after$reducer_fun[[nm]], before$reducer_fun[[nm]], info = nm)
  }
  for (nm in names(before$posthoc_fun)) {
    expect_identical(after$posthoc_fun[[nm]], before$posthoc_fun[[nm]], info = nm)
  }
  expect_setequal(fmrigds::get_reducer("meta:fe")$provides,
                  c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"))
  expect_false(any(grepl("^frg:", fmrigds::list_reducers())))
})

test_that("fmrigds meta:fe still computes the closed-form FE estimate after fmrireg runs (B1)", {
  skip_if_not_installed("fmrigds")
  gd <- .d_demo_gds()
  invisible(fmri_meta(gd, formula = ~ 1, method = "fe", verbose = FALSE))
  set.seed(3)
  P <- 4; S <- 6
  b <- array(rnorm(P * S), c(P, S, 1)); v <- array(runif(P * S, 0.5, 2), c(P, S, 1))
  g <- fmrigds::as_gds(list(beta = b, var = v), space = NULL)
  out <- fmrigds::compute(fmrigds::reduce(fmrigds::as_plan(g), method = "meta:fe"))
  est <- as.numeric(fmrigds::assay(out, "beta_g"))
  w <- 1 / v[, , 1]
  expect_equal(est, rowSums(w * b[, , 1]) / rowSums(w), tolerance = 1e-10)
})

test_that("fmrireg's spatial-FDR posthoc is available under a namespaced name", {
  skip_if_not_installed("fmrigds")
  set.seed(42)
  P <- 50; K <- 2
  z_mat <- matrix(rnorm(P * K), nrow = P, ncol = K)
  grp <- rep(1:5, each = 10)
  g <- fmrigds::as_gds(list(z = array(z_mat, dim = c(P, 1, K))), space = NULL)
  pl <- fmrigds::posthoc(fmrigds::as_plan(g), method = "fmrireg:fdr:spatial",
                         options = list(group = grp))
  q <- fmrigds::assay(fmrigds::compute(pl), "q")
  for (k in 1:K) {
    expect_equal(as.numeric(q[, 1, k]), as.numeric(spatial_fdr(z = z_mat[, k], group = grp)$q))
  }
})
