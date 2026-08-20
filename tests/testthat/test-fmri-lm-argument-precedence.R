# Precedence and conflict semantics for fmri_lm()'s configuration arguments.
#
# Nothing asserted this before, which is why three defects survived: robust_psi
# was unreachable, an explicit robust = FALSE could be overridden, and cfg=
# silently discarded every other configuration argument. Phase 1 of
# docs/plans/fmri-lm-interface-redesign.md.

.robust_norm <- function(...) fmrireg:::.fmri_lm_normalize_robust_options(...)
.ar_norm     <- function(...) fmrireg:::.fmri_lm_normalize_ar_options(...)

# ------------------------------------------------------------------ robust ----

test_that("robust_psi is reachable and selects the psi function", {
  # Regression: `robust` used to default to FALSE rather than NULL, so
  # robust_options$type was always occupied before robust_psi was consulted and
  # the argument had no observable effect anywhere.
  expect_equal(.robust_norm(robust = NULL, robust_psi = "bisquare")$type, "bisquare")
  expect_equal(.robust_norm(robust = NULL, robust_psi = "huber")$type, "huber")
})

test_that("unspecified robust settings default to no robust fitting", {
  expect_false(.robust_norm(robust = NULL)$type)
  expect_false(.robust_norm(robust = FALSE)$type)
})

test_that("legacy robust settings must be scalar and non-missing", {
  expect_error(.robust_norm(robust = c(TRUE, FALSE)), "single non-missing")
  expect_error(.robust_norm(robust = NA), "single non-missing")
  expect_error(.robust_norm(robust = 1), "must be TRUE")
  expect_error(.robust_norm(robust = NULL, robust_psi = c("huber", "bisquare")),
               "single non-missing")
  expect_error(.robust_norm(robust = NULL,
                            robust_options = list(type = c("huber", "bisquare"))),
               "single non-missing")
})

test_that("an explicit robust = FALSE is not silently overridden", {
  # Regression: robust = FALSE + robust_options = list(type = "huber") used to
  # return "huber", i.e. explicitly disabling robust fitting enabled it.
  expect_error(
    .robust_norm(robust = FALSE, robust_options = list(type = "huber")),
    "Conflicting robust settings"
  )
})

test_that("robust settings that agree are accepted", {
  expect_equal(.robust_norm(robust = "huber",
                            robust_options = list(type = "huber"))$type, "huber")
  # TRUE canonicalises to huber, so this agrees rather than conflicts
  expect_equal(.robust_norm(robust = TRUE,
                            robust_options = list(type = "huber"))$type, "huber")
})

test_that("non-type robust options coexist with a robust shorthand", {
  # the shapes that actually occur in this repo's call sites
  out <- .robust_norm(robust = "bisquare", robust_options = list(c_tukey = 4.685))
  expect_equal(out$type, "bisquare")
  expect_equal(out$c_tukey, 4.685)

  out2 <- .robust_norm(robust = FALSE, robust_options = list(max_iter = 10L))
  expect_false(out2$type)
  expect_equal(out2$max_iter, 10L)
})

test_that("disagreeing robust spellings error whichever pair is used", {
  expect_error(.robust_norm(robust = "huber", robust_psi = "bisquare"),
               "Conflicting robust settings")
  expect_error(.robust_norm(robust = NULL, robust_psi = "huber",
                            robust_options = list(type = "bisquare")),
               "Conflicting robust settings")
})

# ---------------------------------------------------------------------- AR ----

test_that("AR shorthands apply when the options list does not set the key", {
  expect_equal(.ar_norm(cor_struct = "ar2")$struct, "ar2")
  expect_equal(.ar_norm(ar_options = list(struct = "arp"), ar_p = 4)$p, 4)
  expect_true(.ar_norm(ar_options = list(struct = "ar1"), ar_voxelwise = TRUE)$voxelwise)
})

test_that("an AR shorthand that disagrees with ar_options is an error", {
  # Regression: the list used to win silently. `ar_voxelwise` was additionally
  # documented as overriding ar_options$voxelwise, which was backwards.
  expect_error(.ar_norm(ar_options = list(struct = "ar2"), cor_struct = "ar1"),
               "Conflicting AR settings")
  expect_error(.ar_norm(ar_options = list(struct = "ar1", voxelwise = FALSE),
                        ar_voxelwise = TRUE),
               "Conflicting AR settings")
  expect_error(.ar_norm(ar_options = list(struct = "arp", p = 2), ar_p = 4),
               "Conflicting AR settings")
})

test_that("an AR shorthand that agrees with ar_options is accepted", {
  expect_equal(.ar_norm(ar_options = list(struct = "ar2"), cor_struct = "ar2")$struct, "ar2")
})

test_that("voxelwise defaults to FALSE when nothing requests it", {
  expect_false(.ar_norm(ar_options = list(struct = "ar1"))$voxelwise)
  expect_false(.ar_norm()$voxelwise)
})

# --------------------------------------------------------------------- cfg ----

test_that("cfg cannot be combined with other configuration arguments", {
  # Regression: .fmri_lm_build_config computed the full normalisation and then
  # discarded it when cfg was present, so cor_struct = "ar2" alongside a cfg
  # naming "iid" silently ran iid.
  cfg <- fmri_lm_control(ar_options = list(struct = "iid"))
  expect_error(
    fmrireg:::.fmri_lm_build_config(robust = NULL, cor_struct = "ar2", engine_cfg = cfg),
    "cannot be combined"
  )
  expect_error(
    fmrireg:::.fmri_lm_build_config(robust = TRUE, engine_cfg = cfg),
    "cannot be combined"
  )
})

test_that("cfg alone is passed through unchanged", {
  cfg <- fmri_lm_control(ar_options = list(struct = "ar2"))
  out <- fmrireg:::.fmri_lm_build_config(robust = NULL, engine_cfg = cfg)
  expect_equal(out$ar$struct, "ar2")
  expect_s3_class(out, "fmri_lm_control")
})

test_that("malformed cfg is rejected instead of silently ignored", {
  expect_error(
    fmrireg:::.fmri_lm_build_config(robust = NULL, engine_cfg = 42),
    "created by `fmri_lm_control\\(\\)`"
  )
  expect_error(
    fmrireg:::.fmri_lm_build_config(robust = NULL, engine_cfg = list()),
    "created by `fmri_lm_control\\(\\)`"
  )
})

# ------------------------------------------------------------------ engine ----

test_that("engine warns about the execution arguments it ignores", {
  w <- fmrireg:::.fmri_lm_warn_engine_ignores
  cl <- quote(fmri_lm(f, dataset = d, engine = "rrr_gls", strategy = "chunkwise",
                      nchunks = 4))
  expect_warning(w(cl, "rrr_gls"), "ignores")
  expect_warning(w(cl, "rrr_gls"), "strategy")
  expect_warning(w(cl, "rrr_gls"), "nchunks")

  # nothing ignorable supplied -> silent
  expect_silent(w(quote(fmri_lm(f, dataset = d, engine = "rrr_gls")), "rrr_gls"))
})

# ------------------------------------------------------- end-to-end on fits ----

test_that("fmri_lm honours each configuration route and rejects conflicts", {
  skip_if_not_installed("fmridataset")
  set.seed(451)
  TR <- 2; Tr <- 100L; n_run <- 2L; n_t <- Tr * n_run
  onsets <- sort(runif(24, 0, Tr * TR - 20))
  cond <- factor(rep(c("a", "b"), length.out = 24))
  etab <- data.frame(onset = rep(onsets, n_run), cond = rep(cond, n_run),
                     run = rep(seq_len(n_run), each = 24))
  dset <- matrix_dataset(matrix(rnorm(n_t * 4), n_t, 4), TR = TR,
                         run_length = rep(Tr, n_run), event_table = etab)
  f <- onset ~ hrf(cond)

  expect_s3_class(fmri_lm(f, block = ~ run, dataset = dset), "fmri_lm")
  expect_s3_class(fmri_lm(f, block = ~ run, dataset = dset, cor_struct = "ar2"), "fmri_lm")
  expect_s3_class(fmri_lm(f, block = ~ run, dataset = dset,
                          ar_options = list(struct = "ar1")), "fmri_lm")
  expect_s3_class(fmri_lm(f, block = ~ run, dataset = dset, robust = FALSE), "fmri_lm")
  expect_s3_class(fmri_lm(f, block = ~ run, dataset = dset,
                          cfg = fmri_lm_control(ar_options = list(struct = "iid"))),
                  "fmri_lm")

  expect_error(
    fmri_lm(f, block = ~ run, dataset = dset, cor_struct = "ar2",
            cfg = fmri_lm_control(ar_options = list(struct = "iid"))),
    "cannot be combined"
  )
  expect_error(
    fmri_lm(f, block = ~ run, dataset = dset, cor_struct = "ar1",
            ar_options = list(struct = "ar2")),
    "Conflicting AR settings"
  )
  expect_error(
    fmri_lm(f, block = ~ run, dataset = dset, robust = FALSE,
            robust_options = list(type = "huber")),
    "Conflicting robust settings"
  )
})
