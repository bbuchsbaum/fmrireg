# Regression tests for the 2026-10 simplification, stream C (engines, plugin
# API, latent fits, betas).

# --- B7: rrr_gls resolves censoring like every other fitting path ---------

# Two runs whose censored volumes hold gross artefacts, so the AR estimate
# moves visibly when censoring is (or is not) applied.
.rrr_censor_fixture <- function(bad = c(10L, 11L, 12L, 55L, 56L)) {
  set.seed(42)
  runs <- c(40L, 40L)
  total <- sum(runs)
  datamat <- matrix(rnorm(total * 6L), total, 6L)
  datamat[bad, ] <- datamat[bad, ] + 50
  ev <- do.call(rbind, lapply(seq_along(runs), function(r) {
    on <- seq(4, runs[r] * 2 - 12, by = 8)
    data.frame(onset = on,
               condition = factor(rep(c("A", "B"), length.out = length(on))),
               run = r)
  }))
  cen <- rep(FALSE, total)
  cen[bad] <- TRUE
  list(
    censored = matrix_frame(datamat, TR = 2, run_length = runs,
                            event_table = ev, censor = cen),
    plain = matrix_frame(datamat, TR = 2, run_length = runs, event_table = ev),
    bad = bad
  )
}

.rrr_censor_fit <- function(ds, censor = NULL) {
  fmri_lm(
    onset ~ hrf(condition), block = ~ run, dataset = ds,
    engine = "rrr_gls", engine_args = list(rank = 1L),
    control = fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec("ar1", censor = censor)
    )
  )
}

.rrr_phi <- function(fit) unlist(fit[["ar_coef"]], use.names = FALSE)
.rrr_coef <- function(fit) as.numeric(fit[["betas_fixed"]])

test_that("rrr_gls uses a frame-carried censor column (B7)", {
  fx <- .rrr_censor_fixture()

  from_column <- .rrr_censor_fit(fx$censored)
  explicit <- .rrr_censor_fit(fx$plain, censor = fx$bad)
  auto <- .rrr_censor_fit(fx$censored, censor = "auto")
  uncensored <- .rrr_censor_fit(fx$plain)

  # The censor column is honoured exactly as the explicit indices are.
  expect_equal(.rrr_phi(from_column), .rrr_phi(explicit), tolerance = 1e-12)
  expect_equal(.rrr_coef(from_column), .rrr_coef(explicit), tolerance = 1e-12)
  expect_equal(.rrr_phi(auto), .rrr_phi(explicit), tolerance = 1e-12)
  expect_equal(.rrr_coef(auto), .rrr_coef(explicit), tolerance = 1e-12)
  # ... and censoring the artefacts changes the AR estimate.
  expect_false(isTRUE(all.equal(.rrr_phi(explicit), .rrr_phi(uncensored))))
})

test_that("rrr_gls censor = 'none' ignores the frame column (B7)", {
  fx <- .rrr_censor_fixture()

  # Before the fix the raw string reached fmriAR as as.integer("none").
  expect_no_warning(
    opted_out <- .rrr_censor_fit(fx$censored, censor = "none")
  )
  uncensored <- .rrr_censor_fit(fx$plain)

  expect_equal(.rrr_phi(opted_out), .rrr_phi(uncensored), tolerance = 1e-12)
  expect_equal(.rrr_coef(opted_out), .rrr_coef(uncensored), tolerance = 1e-12)
})
