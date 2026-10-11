# Performance pass 2, stream P3-1: the meta kernels compute the fixed-effects
# fit once per feature and reuse it for Q / I2, the DerSimonian-Laird moment
# estimator and the Paule-Mandel starting value; DL uses a vectorised leverage.

p3_meta_sim <- function(S = 18L, P = 12L, seed = 7) {
  set.seed(seed)
  covars <- data.frame(
    x = stats::rnorm(S),
    g = factor(rep(c("a", "b", "c"), length.out = S))
  )
  V <- matrix(stats::runif(S * P, 0.02, 0.4), S, P)
  tau <- rep(c(0, 0.2, 0.5, 0.8), length.out = P)
  Y <- matrix(NA_real_, S, P)
  for (j in seq_len(P)) {
    Y[, j] <- 0.2 + 0.4 * covars$x + stats::rnorm(S, sd = sqrt(V[, j] + tau[j]^2))
  }
  # missing subjects: NA effects and NA / non-positive variances
  Y[c(2, 9), 3] <- NA
  V[5, 4] <- NA
  V[c(1, 7, 11), 6] <- 0
  list(Y = Y, V = V, covars = covars)
}

# Independent R reference for the closed-form estimators (FE and DL), and the
# same bisection Paule-Mandel the kernel documents.
p3_ref_fit <- function(y, v, X, method) {
  keep <- is.finite(y) & is.finite(v) & v > 0
  y <- y[keep]; v <- v[keep]; X <- X[keep, , drop = FALSE]
  m <- length(y); p <- ncol(X); df <- m - p
  wls <- function(w) {
    inv <- solve(crossprod(X * sqrt(w)))
    b <- inv %*% crossprod(X, w * y)
    r <- y - X %*% b
    list(beta = drop(b), inv = inv, Q = sum(w * r^2))
  }
  fe <- wls(1 / v)
  tau2 <- 0
  if (method == "dl") {
    w0 <- 1 / v
    lev <- rowSums((X %*% fe$inv) * X)
    C <- sum(w0) - sum(w0^2 * lev)
    tau2 <- max(0, (fe$Q - df) / C)
  } else if (method == "pm" && fe$Q - df > 0) {
    f <- function(t2) wls(1 / (v + t2))$Q - df
    hi <- max(1, 1000 * stats::median(v))
    for (k in 1:20) { if (f(hi) <= 0) break; hi <- hi * 2 }
    a <- 0; b <- hi
    for (it in 1:50) {
      mid <- (a + b) / 2; fm <- f(mid)
      if (abs(fm) < 1e-6) break
      if (fm > 0) a <- mid else b <- mid
    }
    tau2 <- mid
  }
  fin <- wls(1 / (v + tau2))
  list(beta = unname(fin$beta), se = unname(sqrt(diag(fin$inv))), tau2 = tau2, Q = fe$Q, df = df,
       keep = keep)
}

test_that("P3-1: kernels match metafor::rma (FE, DL, PM; intercept and regression; missing subjects)", {
  skip_if_not_installed("metafor")
  sim <- p3_meta_sim()
  designs <- list(
    intercept = list(X = matrix(1, nrow(sim$Y), 1), mods = NULL),
    regression = list(X = stats::model.matrix(~ x + g, sim$covars), mods = ~ x + g)
  )
  methods <- c(fe = "FE", dl = "DL", pm = "PM")
  ctrl <- list(tol = 1e-13, threshold = 1e-13, maxiter = 10000)
  for (dn in names(designs)) {
    d <- designs[[dn]]
    for (m in names(methods)) {
      fit <- meta_fit_cov_cpp(sim$Y, sim$V, d$X, m, "none", 1.345, 2L, 1L)
      for (j in seq_len(ncol(sim$Y))) {
        keep <- is.finite(sim$Y[, j]) & is.finite(sim$V[, j]) & sim$V[, j] > 0
        cv <- sim$covars[keep, , drop = FALSE]
        ref <- if (is.null(d$mods)) {
          metafor::rma(yi = sim$Y[keep, j], vi = sim$V[keep, j], method = methods[[m]],
                       control = ctrl)
        } else {
          metafor::rma(yi = sim$Y[keep, j], vi = sim$V[keep, j], mods = d$mods,
                       data = cv, method = methods[[m]], control = ctrl)
        }
        info <- paste(dn, m, "voxel", j)
        expect_equal(as.numeric(fit$beta[, j]), as.numeric(ref$beta), tolerance = 1e-6, info = info)
        expect_equal(as.numeric(fit$se[, j]), as.numeric(ref$se), tolerance = 1e-6, info = info)
        expect_equal(fit$tau2[j], ref$tau2, tolerance = 1e-6, info = info)
        expect_equal(fit$Q_fe[j], ref$QE, tolerance = 1e-6, info = info)
        expect_equal(fit$df[j], ref$k - ref$p, info = info)
        # packed covariance (row-major upper triangle) equals vcov(rma)
        vc <- stats::vcov(ref)
        K <- ncol(d$X)
        packed <- unlist(lapply(seq_len(K), function(a) vc[a, a:K]))
        expect_equal(as.numeric(fit$cov_tri[, j]), as.numeric(packed), tolerance = 1e-6, info = info)
      }
    }
  }
})

test_that("P3-1: kernels match the R reference of the documented algorithm to 1e-10", {
  sim <- p3_meta_sim(S = 24L, P = 20L, seed = 3)
  Xs <- list(matrix(1, 24, 1), stats::model.matrix(~ x + g, sim$covars))
  for (X in Xs) for (m in c("fe", "dl", "pm")) {
    fit <- meta_fit_cpp(sim$Y, sim$V, X, m, "none", 1.345, 2L, 1L)
    for (j in seq_len(ncol(sim$Y))) {
      ref <- p3_ref_fit(sim$Y[, j], sim$V[, j], X, m)
      info <- paste(m, ncol(X), j)
      expect_equal(as.numeric(fit$beta[, j]), ref$beta, tolerance = 1e-10, info = info)
      expect_equal(as.numeric(fit$se[, j]), ref$se, tolerance = 1e-10, info = info)
      expect_equal(fit$tau2[j], ref$tau2, tolerance = 1e-10, info = info)
      expect_equal(fit$Q_fe[j], ref$Q, tolerance = 1e-10, info = info)
    }
  }
})

test_that("P3-1 guard: exactly one fixed-effects WLS solve per feature", {
  sim <- p3_meta_sim()
  X <- stats::model.matrix(~ x + g, sim$covars)
  for (m in c("fe", "dl", "pm")) {
    for (rb in c("none", "huber")) {
      d <- meta_fit_diag_cpp(sim$Y, sim$V, X, m, rb, 1.345, 2L)
      expect_equal(as.integer(d$n_fe_solves), rep(1L, ncol(sim$Y)), info = paste(m, rb))
    }
  }
  # FE without robustness needs no further solve at all; DL needs none either
  d_fe <- meta_fit_diag_cpp(sim$Y, sim$V, X, "fe", "none", 1.345, 2L)
  expect_equal(as.integer(d_fe$n_solves), rep(1L, ncol(sim$Y)))
  d_dl <- meta_fit_diag_cpp(sim$Y, sim$V, X, "dl", "none", 1.345, 2L)
  dl_fit <- meta_fit_cpp(sim$Y, sim$V, X, "dl", "none", 1.345, 2L, 1L)
  # one FE solve, plus one final solve only when tau2 > 0
  expect_equal(as.integer(d_dl$n_solves), 1L + as.integer(dl_fit$tau2 > 0))
})
