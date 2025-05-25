skip_if_not_installed("bidser")

BF <- fmrireg:::build_fmri_model_from_dsl

# Simple helper to create dummy models
create_dummy_models <- function() {
  sframe <- sampling_frame(blocklens = c(2), TR = 2)
  ev_mod <- event_model(~ hrf(cond),
                        data = data.frame(onset = c(1),
                                           duration = c(1),
                                           cond = factor("a"),
                                           block = 1),
                        block = ~ block,
                        sampling_frame = sframe,
                        durations = 1)
  bl_mod <- baseline_model(basis = "constant", sframe = sframe)
  list(ev_mod = ev_mod, bl_mod = bl_mod)
}

test_that("fmri_model assembled from components", {
  mods <- create_dummy_models()
  fm <- BF(mods$ev_mod, mods$bl_mod)
  expect_s3_class(fm, "fmri_model")
})
