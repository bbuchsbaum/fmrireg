# Issue #263: fmri_latent_lm() silently dropped `control =` (swallowed by
# `...`) and always fitted joint estimation with one AR model shared across
# components.

.latent263_data <- function(seed = 1L) {
  skip_if_not_installed("neuroim2")
  set.seed(seed)
  rl <- c(120L, 120L); n <- sum(rl); k <- 5L; nvox <- 50L
  ev <- data.frame(run = rep(1:2, each = 10),
                   onset = rep(seq(10, 200, length.out = 10), 2),
                   cond = factor(rep(c("a", "b"), 10)))
  mask <- neuroim2::LogicalNeuroVol(
    array(c(rep(TRUE, nvox), rep(FALSE, 14)), c(4, 4, 4)),
    neuroim2::NeuroSpace(c(4, 4, 4))
  )
  # components with very different temporal autocorrelation
  scores <- sapply(c(.05, .3, .6, .8, .9),
                   function(phi) as.numeric(stats::arima.sim(list(ar = phi), n)))
  scores[, 1] <- scores[, 1] + 2 * rep(c(0, 1), length.out = n)
  L <- matrix(rnorm(nvox * k), nvox, k)
  list(
    latent = latent_list_frame(list(basis = scores, loadings = L, mask = mask),
                          TR = 2, run_length = rl, events = ev),
    scores = as_fmri_frame(scores, TR = 2, run_length = rl, events = ev),
    L = L
  )
}

.latent263_ctrl <- function() {
  fmri_lm_control(
    estimation = estimation_spec("runwise_meta"),
    noise = noise_spec("ar1", pooling = "run", voxelwise = TRUE)
  )
}

test_that("#263: control= is honoured and matches fmri_lm on the score matrix", {
  d <- .latent263_data()
  f <- onset ~ hrf(cond)
  ctrl <- .latent263_ctrl()

  fit_default <- fmri_latent_lm(f, block = ~run, dataset = d$latent,
                                autocor = "ar1", durations = 0)
  fit_ctrl <- fmri_latent_lm(f, block = ~run, dataset = d$latent,
                             durations = 0, control = ctrl)
  # autocor that agrees with control is accepted
  fit_both <- fmri_latent_lm(f, block = ~run, dataset = d$latent,
                             durations = 0, autocor = "ar1", control = ctrl)
  ref <- fmri_lm(f, block = ~run, dataset = d$scores, control = ctrl)

  expect_s3_class(fit_ctrl, "fmri_latent_lm")
  # Previously identical: control was discarded.
  expect_gt(max(abs(as.matrix(coef(fit_default)) - as.matrix(coef(fit_ctrl)))), 1e-3)
  expect_equal(as.matrix(coef(fit_ctrl)), as.matrix(coef(ref)), tolerance = 1e-8)
  expect_equal(as.matrix(standard_error(fit_ctrl)), as.matrix(standard_error(ref)),
               tolerance = 1e-8)
  expect_equal(as.matrix(coef(fit_both)), as.matrix(coef(fit_ctrl)))
})

test_that("#263: unknown arguments and conflicting autocor are errors", {
  d <- .latent263_data()
  f <- onset ~ hrf(cond)
  expect_error(
    fmri_latent_lm(f, block = ~run, dataset = d$latent, contorl = .latent263_ctrl()),
    "Unused argument.*contorl"
  )
  expect_error(
    fmri_latent_lm(f, block = ~run, dataset = d$latent, autocor = "ar2",
                   control = .latent263_ctrl()),
    "conflicts with the noise model"
  )
  expect_error(
    fmri_latent_lm(f, block = ~run, dataset = d$latent, control = list()),
    "fmri_lm_control"
  )
  expect_error(
    fmri_latent_lm(f, block = ~run, dataset = d$latent,
                   control = fmri_lm_control(robust = robust_spec("huber"))),
    "Robust fitting is not supported"
  )
})

test_that("#263: durations defaults to 0 like fmri_lm()", {
  d <- .latent263_data()
  f <- onset ~ hrf(cond)
  a <- fmri_latent_lm(f, block = ~run, dataset = d$latent)
  b <- fmri_latent_lm(f, block = ~run, dataset = d$latent, durations = 0)
  expect_equal(as.matrix(coef(a)), as.matrix(coef(b)))
})

test_that("#263: per-component fits refuse voxel-space SEs but reconstruct betas", {
  d <- .latent263_data()
  f <- onset ~ hrf(cond)
  fit <- fmri_latent_lm(f, block = ~run, dataset = d$latent,
                        control = .latent263_ctrl())
  B <- as.matrix(coef(fit))
  expect_equal(unname(as.matrix(coef(fit, recon = TRUE))), unname(d$L %*% B),
               tolerance = 1e-10)
  expect_error(standard_error(fit, recon = TRUE), "shared across components")
  expect_error(stats(fit, recon = TRUE), "shared across components")

  # The default shared-noise fit still supports voxel-space inference.
  fit_joint <- fmri_latent_lm(f, block = ~run, dataset = d$latent, autocor = "ar1")
  se <- standard_error(fit_joint, recon = TRUE)
  expect_equal(nrow(se), nrow(d$L))
  expect_true(all(is.finite(as.matrix(se))))
})
