devtools::load_all()

dset <- readRDS("data-raw/hyper_roi_1006.rds")

mats <- dset$mats[1:3]
conf <- dset$conf[1:3]

xmats <- lapply(1:length(mats), function(i) {
  X <- mats[[i]]
  nuisance <- conf[[i]]
  nuisance <- apply(nuisance, 2, function(vals) {
    vals[is.na(vals)] <- 0
    vals
  })
  
  nuisance <- scale(nuisance)
  lm.1 <- lm(X ~ nuisance[,c(1,2,4,8,9,10, 20:23)])
  Xr <- resid(lm.1)
})

evtab <- do.call(rbind, dset$des[1:3])
evtab$onset <- evtab$Onset/1000
evtab$stim <- factor(rep(paste0("stim", 1:6), length.out=nrow(evtab)))
mdset <- matrix_dataset(do.call(rbind, xmats), 1.77, rep(210,3), evtab)


params <- expand.grid(h1=seq(1,3, by=1),  h2=seq(3,7,by=1), h3=seq(6,8,by=1), h4=seq(5,9,by=1), f1=seq(0,.2,by=.1), f2=seq(0,.2,by=.1))
p=gen_hrf_library(hrf_half_cosine, params)
pmat <- p(0:24)
L <- pmat
pres <- multivarious::pca(pmat, preproc=multivarious::pass())
basis <- pres$u[,1:6]
Lk <- basis

basis_set <- lapply(1:ncol(basis), function(i) { gen_empirical_hrf(0:24, basis[,i])})
hrf_list <- do.call(gen_hrf_set, basis_set)


evtab$stim <- factor(rep(paste0("stim", 1:4, length.out=nrow(evtab))))
levs <- levels(evtab$stim)
mdset <- matrix_dataset(do.call(rbind, xmats), 1.77, rep(210,3), evtab)
con1 <- unit_contrast(~ stim %in% levs, where= ~ basis == "basis01", name="stim_baseline")
con1 <- unit_contrast(~ stim %in% levs, name="stim_baseline")
emphrf <- gen_empirical_hrf(0:24, pmat[,1576])
emod <- event_model(onset ~ hrf(stim, basis=hrf_list, contrasts=con1), data=mdset$event_table, block = ~ run, 
                    sampling_frame=mdset$sampling_frame)
emod <- event_model(onset ~ hrf(stim, basis=emphrf, contrasts=con1), data=mdset$event_table, block = ~ run, 
                    sampling_frame=mdset$sampling_frame)
bmod <- baseline_model(basis="constant", sframe=mdset$sampling_frame)
dmat <-design_matrix(emod)
fmod <- fmri_model(emod, bmod)
fit <- fmri_lm_fit(fmod, mdset)
con = stats(fit, "contrasts")
betas <- coef(fit)


dmat <- design_matrix(emod)

bmat <- as.matrix(dmat)

X <- scale(mdset$datamat)
pres4 <- multivarious::pca(X, preproc=multivarious::pass())
Mk <- pres4$u[,1:100]
Rk <- design_matrix(emod)


### key computations
C <- t(Rk) %*% Mk

## average covarance matrix across stimulus groups
C <- multivarious::group_means(Y=rep(1:6, 4), C)
svd_result <- svd(C)

Xproj <- t(X) %*% Mk %*% svd_result$v
Lproj <- t(L) %*% Lk %*% svd_result$u

basis_euc <- proxy::simil(Lproj, Xproj, method="euclidean")
best = apply(basis_euc, 2, function(x) which.max(abs(x)))
####


C <- cov(bmat, smat)
svd.1 <- svd(t(smat) %*% bmat)
Xreduced <- svd.1$u[1:6,] %*% t(smat)
basis_proj <- t(pmat) %*% basis %*% svd.1$u
Xproj <- smat %*% svd.1$v

basis_euc <- proxy::simil(basis_proj, svd.1$v, method="euclidean")
best = apply(basis_euc, 2, function(x) which.max(abs(x)))






ons <- global_onsets.sampling_frame(mdset$sampling_frame, evtab$onset,evtab$run)
cvals <- unlist(lapply(1:ncol(pmat), function(i) {
  reg1 <- regressor(onsets=ons, hrf=gen_empirical_hrf(0:24, pmat[,i]))
  x1 <- evaluate(reg1, samples(mdset$sampling_frame, global=TRUE))
  cor(x1, mdset$datamat[,328])
}))


regressors <- lapply(1:ncol(pmat), function(i) {
  print(i)
  reg1 <- regressor(onsets=ons, hrf=gen_empirical_hrf(0:24, pmat[,i]))
  x1 <- evaluate(reg1, samples(mdset$sampling_frame, global=TRUE))
  x1
})

tmp <- lapply(1:ncol(mdset$datamat), function(j) {
  print(j)
  cvals <- unlist(lapply(1:ncol(pmat), function(i) {
    #reg1 <- regressor(onsets=ons, hrf=gen_empirical_hrf(0:24, pmat[,i]))
    #x1 <- evaluate(reg1, samples(mdset$sampling_frame, global=TRUE))
    cor(regressors[[i]], mdset$datamat[,j])
  }))
})
tmp = do.call(cbind, tmp)



bmat <- design_matrix(bmod)

y <- mdset$datamat[,1]
lm.1 <- lm(y ~ as.matrix(dmat2) + as.matrix(bmat)-1)
out <- lapply(1:ncol(mdset$datamat), function(i) {
  print(i)
  y <- mdset$datamat[,i]
  lm.1 <- lm(y ~ as.matrix(dmat2)-1)
  #anova(lm.1)$F[1:2]
  coef(summary(lm.1))[1:6,3]
})
out2 <- do.call(rbind, out)

bmat <- as.matrix(dmat)
smat <- scale(mdset$datamat)
C <- cov(bmat, smat)
svd.1 <- svd(t(smat) %*% bmat)
Xreduced <- svd.1$u[1:6,] %*% t(smat)
basis_proj <- t(pmat) %*% basis %*% svd.1$u
Xproj <- smat %*% svd.1$v

basis_euc <- proxy::simil(basis_proj, svd.1$v, method="euclidean")
best = apply(basis_euc, 2, function(x) which.max(abs(x)))


L <- pmat
pres1 <-  multivarious::pca(L, preproc=multivarious::pass())
Lk <- pres$u[,1:6]
basis_set <- lapply(1:ncol(basis), function(i) { gen_empirical_hrf(0:24, Lk[,i])})
hrf_list <- do.call(gen_hrf_set, basis_set)
emod <- event_model(onset ~ hrf(stim, basis=hrf_list), data=mdset$event_table, block = ~ run, sampling_frame=mdset$sampling_frame)
Rk <- design_matrix(emod)
Rk <- apply(Rk, 2, function(vals) vals/sqrt(sum(vals^2)))
X <- mdset$datamat
pres2 <- multivarious::pca(X, preproc=multivarious::standardize())
Xk <- pres2$u[,1:6]

C <- cbind(Rk, Xk)
pres3 <- multivarious::pca(C, preproc=multivarious::pass())

Lproj <- multivarious::partial_project(pres2, Lk, 1:6)
Mproj <- multivarious::partial_project(pres2, Xk, 7:12)

S <- Lproj %*% t(Mproj)


pres4 <- multivarious::pca(X, preproc=multivarious::standardize())
Mk <- pres4$u[,1:67]
C <- t(Rk) %*% Mk
svd_result <- svd(C)

Xproj <- t(X) %*% Mk %*% svd_result$v
Lproj <- t(L) %*% Lk %*% svd_result$u

basis_euc <- proxy::simil(Lproj, Xproj, method="cosine")
best = apply(basis_euc, 2, function(x) which.max(abs(x)))



U <- svd_result$u
S <- svd_result$d
S_inv <- diag(1 / diag(S))
z <- t(Mk) %*% X

C <- Mk %*% t(Rk)
svd_result <- svd(C)

Lproj <- t(Lk) %*% svd_result$v[,1:6]
Xproj <- Mk %*% t(svd_result$u[,1:6])

## L = 25 by 2025
## Lk = 25 by 6
## Rk = 630 by 6
## v = 630 by 6

## t(L) %*% Lk = 2025 by 6




