# Regression tests for simplification stream E (2026-10): print methods,
# with_package(), BIDS export, imports.

.e_two_run_frame <- function() {
  set.seed(11)
  n <- 60
  ev <- data.frame(
    onset = rep(c(5, 20, 35, 50), 2),
    condition = factor(rep(c("A", "B", "C", "A"), 2)),
    run = rep(1:2, each = 4)
  )
  matrix_frame(matrix(rnorm(2 * n * 6), 2 * n, 6), TR = 1,
               run_length = c(n, n), event_table = ev)
}

test_that("print.fmri_lm reports the fit's scope, stored df, formula and contrasts (B2)", {
  dset <- .e_two_run_frame()
  cons <- list(
    fmridesign::pair_contrast(~ condition == "A", ~ condition == "B", name = "AvB"),
    fmridesign::oneway_contrast(~ condition, name = "main_cond")
  )
  fit <- fmri_lm(onset ~ hrf(condition, contrasts = cons), block = ~run,
                 dataset = dset)
  out <- capture.output(print(fit))

  df_stored <- unique(fit[["result"]][["df"]][["inference"]])
  expect_length(df_stored, 1L)
  expect_true(any(grepl("Estimation scope: runwise_meta", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Inference df (residual): %s", format(df_stored)),
                        out, fixed = TRUE)))
  expect_true(any(grepl("onset ~ hrf(condition, contrasts = cons)", out, fixed = TRUE)))
  expect_true(any(grepl("Simple contrasts (1): AvB", out, fixed = TRUE)))
  expect_true(any(grepl("F-contrasts (1): main_cond", out, fixed = TRUE)))
  expect_true(any(grepl("Parameters: 3 event + 8 baseline", out, fixed = TRUE)))

  # joint (chunkwise) scope, AR and robust configuration lines
  fit_ar <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                    strategy = "chunkwise", nchunks = 2, cor_struct = "ar1")
  out_ar <- capture.output(print(fit_ar))
  expect_true(any(grepl("Estimation scope: joint", out_ar, fixed = TRUE)))
  expect_true(any(grepl("AR structure: ar1", out_ar, fixed = TRUE)))
  expect_true(any(grepl(
    sprintf("Inference df (residual): %s", format(unique(fit_ar[["result"]][["df"]][["inference"]]))),
    out_ar, fixed = TRUE)))

  fit_rob <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                     robust = "huber")
  out_rob <- capture.output(print(fit_rob))
  expect_true(any(grepl("Robust method: huber (tuning constant 1.345)", out_rob, fixed = TRUE)))

  # fmri_rlm keeps its banner and falls through to the fmri_lm summary
  fit_rlm <- fmri_rlm(onset ~ hrf(condition), block = ~run, dataset = dset)
  out_rlm <- capture.output(print(fit_rlm))
  expect_true(any(grepl("Robust fMRI Linear Model", out_rlm, fixed = TRUE)))
  expect_true(any(grepl("Estimation scope: runwise_meta", out_rlm, fixed = TRUE)))
})

test_that("print.fmri_lm names the engine and its df for engine fits (B2)", {
  dset <- .e_two_run_frame()
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 engine = "latent_sketch")
  out <- capture.output(print(fit))
  df_stored <- unique(fit[["result"]][["df"]][["inference"]])
  expect_true(any(grepl("Engine: latent_sketch", out, fixed = TRUE)))
  expect_true(any(grepl(
    sprintf("Inference df (%s): %s", fit[["result"]][["df"]][["method"]],
            format(df_stored, digits = 5L)),
    out, fixed = TRUE)))
})

test_that("print.fmri_model shows event terms, event count, baseline basis and total columns (B3)", {
  dset <- .e_two_run_frame()
  base <- baseline_model(basis = "poly", degree = 2,
                         sframe = fmrihrf::sampling_frame(c(60, 60), TR = 1))
  fm <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                          baseline_model = base, dataset = dset)
  out <- capture.output(print(fm))

  expect_true(any(grepl("Formula: onset ~ hrf(condition)", out, fixed = TRUE)))
  expect_true(any(grepl("condition (<event_term>)", out, fixed = TRUE)))
  expect_true(any(grepl("Number of Events: 8", out, fixed = TRUE)))
  expect_true(any(grepl("Basis type: poly", out, fixed = TRUE)))
  expect_true(any(grepl("Degree: 2", out, fixed = TRUE)))
  # 3 condition columns + (2 poly + 1 intercept) per run x 2 runs
  expect_identical(ncol(design_matrix(fm)), 9L)
  expect_true(any(grepl("Total design columns: 9", out, fixed = TRUE)))
})

test_that("with_package checks every package and names all missing ones (B6)", {
  expect_null(fmrireg:::with_package(c("stats", "utils")))
  err <- tryCatch(
    fmrireg:::with_package(c("stats", "notapkg.zzz1", "notapkg.zzz2")),
    error = function(e) conditionMessage(e)
  )
  expect_identical(
    err,
    "Please install the notapkg.zzz1, notapkg.zzz2 packages to use this functionality"
  )
  # single missing package keeps the original wording
  expect_error(fmrireg:::with_package("notapkg.zzz1"),
               "Please install the notapkg.zzz1 package to use this functionality",
               fixed = TRUE)
})
