test_that("2D steerable set and basis reconstruction via LS calibration", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)

  N <- fmrireg.cca:::steerable_dirs_2d()  # 2×3
  expect_equal(dim(N), c(2L,3L))
  for (i in 1:3) expect_equal(sum(N[,i]^2), 1, tolerance = 1e-12)

  # Impulse at center to match kernel-domain identities
  nx <- 41; ny <- 37
  cx <- ceiling(nx/2); cy <- ceiling(ny/2)
  img <- matrix(0, nx, ny); img[cx, cy] <- 1

  spacing <- c(2,2); fwhm <- 6
  base <- fmrireg.cca:::friman_base_responses_2d(img, spacing, fwhm)
  G0  <- base$G0; Gxx <- base$Gxx; Gyy <- base$Gyy; Gxy <- base$Gxy

  # Calibrate A,B for sigma in voxel units
  c <- 1/(2*sqrt(2*log(2)))
  sigma_vox <- (fwhm * c) / spacing[1]
  ab <- fmrireg.cca:::friman_calibrate_AB_2d(sigma_vox, 1.0, 0.0, -1)
  A <- ab$A; B <- ab$B

  # Oriented sum using calibrated A,B
  D2 <- function(nx,ny) nx*nx*Gxx + ny*ny*Gyy + 2*nx*ny*Gxy
  Ori_sum <- (A*G0*3) + B*(D2(N[1,1],N[2,1]) + D2(N[1,2],N[2,2]) + D2(N[1,3],N[2,3]))

  # Expected oriented sum = (1 - g_iso) * G0 (Eq. 11 summed over 3 dirs)
  sigma_half <- sigma_vox/2
  dx <- matrix(rep(1:nx, times=ny), nx, ny) - cx
  dy <- matrix(rep(1:ny, each=nx), nx, ny) - cy
  giso <- exp(-0.5 * (dx*dx + dy*dy)/(sigma_half*sigma_half))
  Ori_target <- (1 - giso) * G0

  # Check oriented approximation quality (mean abs error relative to mean G0)
  mae <- mean(abs(Ori_sum - Ori_target))
  expect_lt(mae/mean(G0), 2e-3)

  # Full reconstruction: center + oriented ≈ G0
  recon <- giso * G0 + Ori_sum
  expect_lt(mean(abs(recon - G0))/mean(G0), 2e-3)
})

test_that("3D steerable directions are normalized", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  N <- fmrireg.cca:::steerable_dirs_3d()
  expect_equal(dim(N), c(3L,6L))
  for (i in 1:6) expect_equal(sum(N[,i]^2), 1, tolerance = 1e-12)
})

test_that("3D oriented sum matches (1 - g_iso) * G0 via LS calibration", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)
  # 3D impulse at center
  nx <- 25; ny <- 23; nz <- 17
  cx <- ceiling(nx/2); cy <- ceiling(ny/2); cz <- ceiling(nz/2)
  vol <- array(0, c(nx,ny,nz)); vol[cx, cy, cz] <- 1
  spacing <- c(2,2,2); fwhm <- 6
  base <- fmrireg.cca:::friman_base_responses_3d(vol, spacing, fwhm)
  G0  <- base$G0; Gxx <- base$Gxx; Gyy <- base$Gyy; Gzz <- base$Gzz
  Gxy <- base$Gxy; Gxz <- base$Gxz; Gyz <- base$Gyz

  # Calibrate A,B using avg sigma(vox) and an arbitrary direction (isotropic)
  c <- 1/(2*sqrt(2*log(2)))
  sigma_vox <- (fwhm * c) / spacing[1]
  ab <- fmrireg.cca:::friman_calibrate_AB_3d(sigma_vox, c(1,0,0), -1)
  A <- ab$A; B <- ab$B

  N3 <- fmrireg.cca:::steerable_dirs_3d()  # 3×6
  # Sum over 6 directions: Ori_sum ≈ (1 - g_iso) * G0
  D2n <- function(nx,ny,nz) nx*nx*Gxx + ny*ny*Gyy + nz*nz*Gzz + 2*nx*ny*Gxy + 2*nx*nz*Gxz + 2*ny*nz*Gyz
  Ori_sum <- 6*A*G0
  for (i in 1:6) {
    Ori_sum <- Ori_sum + B * D2n(N3[1,i], N3[2,i], N3[3,i])
  }

  # Target oriented component: (1 - g_iso) * G0 with sigma/2 center weight
  sigma_half <- sigma_vox/2
  dx <- array(rep(1:nx, times=ny*nz), c(nx,ny,nz)) - cx
  dy <- array(rep(rep(1:ny, each=nx), times=nz), c(nx,ny,nz)) - cy
  dz <- array(rep(1:nz, each=nx*ny), c(nx,ny,nz)) - cz
  giso <- exp(-0.5 * (dx*dx + dy*dy + dz*dz)/(sigma_half*sigma_half))
  Ori_target <- (1 - giso) * G0

  # Relative MAE thresholds
  mae_oriented <- mean(abs(Ori_sum - Ori_target)) / mean(G0)
  expect_lt(mae_oriented, 4e-3)

  recon <- giso * G0 + Ori_sum
  mae_recon <- mean(abs(recon - G0)) / mean(G0)
  expect_lt(mae_recon, 4e-3)
})
