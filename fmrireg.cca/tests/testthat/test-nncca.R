test_that("exact non-negative CCA equals brute-force in tiny dims", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)
  set.seed(42)
  k <- 3; p <- k; q <- 2
  for (trial in 1:10) {
    A <- matrix(rnorm(1000), 50, p); B <- matrix(rnorm(1000), 50, q)
    Cxx <- crossprod(A)/nrow(A) + diag(1e-3, p)
    Cyy <- crossprod(B)/nrow(B) + diag(1e-3, q)
    Cxy <- crossprod(A,B)/nrow(A)
    sol <- fmrireg.cca:::friman_nncca2xk(Cxx, Cyy, Cxy, shrink = 1e-3,
                                         simplex_x = FALSE, simplex_y = FALSE)
    rho <- sol$rho; wx <- sol$wx; wy <- sol$wy
    expect_true(all(wx >= -1e-10)); expect_true(all(wy >= -1e-10))

    best <- -Inf
    for (maskx in 1:((1<<p)-1)) for (masky in 1:((1<<q)-1)) {
      Ix <- which(as.logical(intToBits(maskx)[1:p]))
      Iy <- which(as.logical(intToBits(masky)[1:q]))
      U <- solve(chol(Cxx[Ix,Ix])) %*% Cxy[Ix,Iy] %*% solve(t(chol(Cyy[Iy,Iy])))
      sv <- svd(U); s <- sv$d[1]
      wx0 <- rep(0,p); wy0 <- rep(0,q)
      wx0[Ix] <- backsolve(t(chol(Cxx[Ix,Ix])), sv$u[,1])
      wy0[Iy] <- backsolve(t(chol(Cyy[Iy,Iy])), sv$v[,1])
      if (all(wx0 >= -1e-10) && all(wy0 >= -1e-10)) best <- max(best, s)
    }
    expect_equal(rho, best, tolerance = 1e-7)
  }
})

