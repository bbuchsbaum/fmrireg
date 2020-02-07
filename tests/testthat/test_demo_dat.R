library(purrr)
library(assertthat)
library(dplyr)

dat <- readRDS(system.file("extdata", "auditory/aud_demo_dat.RDS", package = "fmrireg"))


cvars0 <- names(dat$confounds[[1]])[1:28]

pick_confounds <- function(conf, cvars=cvars0, npcs=6) {
  lapply(conf, function(x) {
    sel <- select(x, cvars)
    cnames <- names(sel)
    selm <- as.matrix(sel)
    colnames(selm) <- cvars
    mat <- apply(selm, 2, function(x) {
      med <- median(x, na.rm=TRUE)
      x[is.na(x)] <- med
      x
    })
    
    pres <- prcomp(scale(mat))$x[,1:npcs]
    
  })
}

test_that("demo1", {
  
  des <- dat$design
  des$ons1 <- des$onset1/1000
  des$ons2 <- des$onset2/1000
  des$ons3 <- des$onset3/1000
  des$ons4 <- des$onset4/1000
  des$constant <- rep(1, nrow(des))
  sframe <- sampling_frame(rep(197,6), TR=1.5)
  isi <- des$onset1[2:nrow(des)] - des$onset4[1:(nrow(des)-1)]
  isi[isi < 0] <- median(isi)
  des$isi <- c(0, scale(isi))
  
  g = expand.grid(lag=seq(.01,5, by=.2), width=seq(0.1,5, by=.2))
  gres <- do.call(cbind, lapply(1:nrow(g), function(i) {
    h <- gen_hrf(hrf_spmg1, lag=g[i,1], width=g[i,2])
    h(0:24)
  }))
  
  pres = prcomp(gres)
  p1=gen_empirical_hrf(0:24, pres$x[,1])
  p2=gen_empirical_hrf(0:24, pres$x[,2])
  p3=gen_empirical_hrf(0:24, pres$x[,3])
  p4=gen_empirical_hrf(0:24, pres$x[,4])
  hset <- gen_hrf_set(p1,p2,p3,p4)
 
  
  for (i in 0:5) {
    nuis <- pick_confounds(dat$confounds, npcs=12)
    bmod <- baseline_model(basis="bs", degree=4, nuisance_list=nuis, sframe=sframe)
    h <- gen_hrf(hrf_gaussian, mean=6,sd=2)
    ##hf <- hrf_blocked(width=5, half_life=i, summate=FALSE)
    
    hrf1 <- gen_hrf(hrf_gaussian, lag=0)
    hrf2 <- gen_hrf(hrf_gaussian, lag=3)
    
    emod <- event_model(ons1 ~ hrf(voicing) + hrf(vowel) + hrf(place) + hrf(speaker) + hrf(gender),
                        #hrf(constant, basis=hrf2, prefix="h4") + 
                        #hrf(Poly(isi,3),basis=hrf1, prefix="h1"), 
                        #hrf(Poly(isi,3),basis=hrf2, prefix="h4"),
                        block = ~ block, data=des, 
                        sampling_frame=sframe)
    
    emod <- event_model(ons1 ~ hrf(constant, basis=hset, prefix="h1"),
                          #hrf(constant, basis=hrf2, prefix="h4") + 
                          #hrf(Poly(isi,3),basis=hrf1, prefix="h1"), 
                          #hrf(Poly(isi,3),basis=hrf2, prefix="h4"),
                      block = ~ block, data=des, 
                      sampling_frame=sframe)
    
    
    emod2 <- event_model(ons1 ~ trialwise("tw", basis=hrf_time),
                        block = ~ block, data=des, 
                        sampling_frame=sframe)
  
    y <- apply(dat$mat, 1, function(x) median(x, na.rm=TRUE))
    dmat_ev <- design_matrix(emod)
    dmat_bv <- design_matrix(bmod)
    evmat <- as.matrix(dmat_ev)
    bvmat <- as.matrix(dmat_bv)
    colnames(evmat) <- names(dmat_ev)
    colnames(bvmat) <- names(dmat_bv)
  
    lm.0 <- lm(y ~ bvmat -1)
    lm.00 <- lm(y ~ bvmat[,1:6] + ccres$scores)
    lm.1 <- lm(y ~ evmat + bvmat -1)
    y0 <- resid(lm.0)
    y00 <- resid(lm.00)
    ares <- auto.arima(y0)
    ares00 <- auto.arima(y00)
    ares2 <- arima(y, order=c(2,0,2), xreg=cbind(evmat, bvmat))
    
    print(cor(evmat[,1],resid(lm.0)))
  }
  
})