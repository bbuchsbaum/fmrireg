test_that("steerable 2D not worse than fixed Gaussian on synthetic", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)
  set.seed(123)
  nx <- 32; ny <- 32; T <- 80
  spacing <- c(2,2); fwhm <- 6
  boxcar <- rep(rep(c(1,0), each=8), length.out=T)
  vols <- vector("list", T); for (t in 1:T) vols[[t]] <- matrix(rnorm(nx*ny, sd=1), nx, ny)
  act <- matrix(FALSE, nx, ny); act[12:22, 10:26] <- TRUE
  for (t in 1:T) if (boxcar[t]==1) vols[[t]][act] <- vols[[t]][act] + 0.6

  # Fixed Gaussian baseline: corr with boxcar
  base <- function(slice) fmrireg.cca:::friman_base_responses_2d(slice, spacing, fwhm)$G0
  corr_fixed <- matrix(NA_real_, nx, ny)
  for (x in 1:nx) for (y in 1:ny) {
    ts <- vapply(vols, function(z) base(z)[x,y], numeric(1))
    corr_fixed[x,y] <- suppressWarnings(stats::cor(ts, boxcar))
  }

  V <- nx*ny; mask_lin <- seq_len(V); mask_z <- rep(1L, V)
  w_dir   <- matrix(c(0.6, 0.3, 0.1), V, 3, byrow = TRUE)
  w_step2 <- matrix(c(0.7, 0.3), V, 2, byrow = TRUE)
  acc <- fmrireg.cca:::friman_pass2_xts_sts_2d(vols, spacing, fwhm, cbind(boxcar), w_dir, w_step2, mask_lin, mask_z)
  XtX <- matrix(sum(boxcar^2), 1, 1)
  beta <- matrix(acc$XtS[1,]/XtX[1,1], nx, ny)

  roc_auc <- function(score, truth) {
    s <- as.numeric(score); y <- as.integer(truth)
    o <- order(-s); y <- y[o]
    tpr <- cumsum(y)/sum(y); fpr <- cumsum(1-y)/sum(1-y)
    # Trapezoidal AUC on stepwise ROC curve
    sum(diff(c(0,fpr)) * c(0, head(tpr, -1)))
  }
  auc_fixed  <- roc_auc(corr_fixed, act)
  auc_steer  <- roc_auc(beta, act)
  expect_gte(auc_steer, auc_fixed - 1e-6)
})

