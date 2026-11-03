# Plot labels come from the model's own variable names, not from splitting
# column names on "_".

test_that("factor names containing underscores are not cut apart", {
  vars <- c("trial_type", "rt", "cue", "load")
  expect_equal(.strip_level_prefixes("trial_type.face", vars), "face")
  expect_equal(.strip_level_prefixes("trial_type.face_rt", vars), paste0("face ", .glyph("times"), " rt"))
  expect_equal(.strip_level_prefixes("cue.a_load.hi", vars), paste0("a ", .glyph("times"), " hi"))
  expect_equal(.strip_level_prefixes("cond.a_load.hi"), paste0("a ", .glyph("times"), " hi"))
})

test_that("a constant drift basis is labelled and grouped as run intercepts", {
  sf <- fmrihrf::sampling_frame(c(40, 40), TR = 2)
  ev <- data.frame(onset = c(10, 30, 10, 30), trial_type = c("face", "house", "house", "face"),
                   run = c(1, 1, 2, 2))
  em <- fmridesign::event_model(onset ~ hrf(trial_type), data = ev, block = ~ run, sampling_frame = sf)
  bm <- fmridesign::baseline_model(basis = "constant", sframe = sf)
  dset <- matrix_frame(matrix(stats::rnorm(80 * 3), 80, 3), TR = 2, run_length = c(40, 40), event_table = ev)
  mod <- fmri_model(em, bm, dset)
  info <- .column_info(mod)
  base <- info$kind != "event"
  expect_true(all(info$kind[base] == "intercept"))
  expect_equal(sub(paste0(" ", .glyph("dot"), " run [0-9]+$"), "", info$label[base]), c("intercept", "intercept"))
  expect_equal(info$label[!base], c("face", "house"))
})

test_that("palettes and type follow the shared neuro.* options", {
  withr::local_options(neuro.palette.categorical = NULL, neuro.palette.diverging = NULL)
  expect_equal(.fmrireg_categorical(2), c("#0072B2", "#D55E00"))
  expect_equal(.fmrireg_diverging()[["high"]], "#B2182B")
  withr::local_options(neuro.palette.categorical = c("#8a6d1f", "#2f7d6d"),
                       neuro.palette.diverging = c("#2a5bd7", "#f6f5f1", "#d33f29"),
                       neuro.base_family = "serif")
  expect_equal(.fmrireg_categorical(3), c("#8a6d1f", "#2f7d6d", "#8a6d1f"))
  expect_equal(unname(.fmrireg_diverging()), c("#2a5bd7", "#f6f5f1", "#d33f29"))
  expect_equal(theme_fmrireg()$text$family, "serif")
})
