# 
# hrfg <- function(time, amp, mean, sd, c) {
#   amp*dnorm(time, mean=mean, sd=sd) + c
# }
# 
# fit_nls <- function(y, onsets, blockids, bmodel, sframe) {
#   bmat <- design_matrix(bmat)
#   lm.0 <- lm(y ~ bmat)
#   yr <- resid(lm.0)
# 
#   c <- mean(y0)
#   efac <- EV(rep(1, length(onsets)), onsets=onsets,  blockids=blockids, name="ev")
#   eterm <- event_term(list(fac=rep(1,length(onsets))), onsets=onsets, blockids=blockids, durations=0)
#   convolve(eterm, hrf_gaussian,sampling_frame=sframe)
#   reg <- regressor(onsets, hrf=hrf_gaussian)
# 
#   lm.2 <- lm()
# 
# }


library(fmrireg)

facedes <- read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), header=TRUE)
facedes$repnum <- factor(facedes$rep_num)
facedes <- subset(facedes, repnum != "-1")
facedes$repnum <- droplevels(facedes$repnum)

sframe <- sampling_frame(rep(430/2,6), TR=2)
ev <- event_model(onset ~ hrf(repnum, basis="gaussian"), data=facedes, block= ~ run, 
                  sampling_frame=sframe)


gonsets <- global_onsets(sframe, facedes$onset, blockids(ev))

hg <- function(t) {
  hrf_gaussian(t, mean=6.7, sd=1.6)
}

#reg0_signal <- regressor(gonsets[facedes$rep_num == -1], hrf=gen_hrf(hg, span=24))
reg1_signal <- regressor(gonsets[facedes$repnum == 1], hrf=gen_hrf(hg, span=24))
reg2_signal <- regressor(gonsets[facedes$repnum == 2], hrf=gen_hrf(hg, span=24))
reg3_signal <- regressor(gonsets[facedes$repnum == 3], hrf=gen_hrf(hg, span=24))
reg4_signal <- regressor(gonsets[facedes$repnum == 4], hrf=gen_hrf(hg, span=24))

gsam <- samples(sframe, global=TRUE)
y1 <- evaluate(reg1_signal, gsam)
y2 <- evaluate(reg2_signal, gsam)
y3 <- evaluate(reg3_signal, gsam)
y4 <- evaluate(reg4_signal, gsam)
y <- y1+y2+y3+y4
y <- y + rnorm(length(y), sd=.02)
                         
fit_nls <- function(y, event_model, hrf=c("gaussian", "half_cos"), span=24, resolution=.25) {
  term <- terms(event_model)[[1]]
  if (is_continuous(term$evterm)) {
    stop("no categorical terms found: fit_nls currently only works for categorical factors.")
  }
  
  hrf <- match.arg(hrf)
  hrf_fun <- if (hrf == "gaussian") {
    bounds <- list(
      lower=c(1,3),
      upper=c(4,9))
    HRF_GAUSSIAN
   
  } else {
    hrfp <- purrr::partial(hrf_half_cosine, h1=.2, f1=.01, f2=0.1)
    bounds <- list(
      lower=c(2,2,2),
      upper=c(9,9,9))
    gen_hrf(hrfp, span=span)
    
  }
 
  reglist <- regressors.event_term(term$evterm, hrf_fun, event_model$sampling_frame)
  gsamples <- samples(event_model$sampling_frame, global=TRUE)
  
  implist <- lapply(reglist, function(reg) {
    fmrireg:::neural_input(reg, gsamples[1], gsamples[length(gsamples)], resolution=resolution)
  })
  
  gen_ith_prediction <- function(i, a, params) {
    h <- gen_hrf(function(t) {
      do.call(hrf_fun, c(list(t), unlist(params)))
    }, span=span)
    
    hsam <- h(seq(0,span, by=resolution))
    input <- as.numeric(implist[[i]]$neural_input)
    
    highres <- pracma::conv(as.numeric(implist[[i]]$neural_input), hsam)
    t <- implist[[i]]$time
    out <- approx(t, highres[1:length(t)], xout=gsamples)$y
    #reg <- regressor(reglist[[i]]$onsets, hrf=h, amplitude=rep(a, length(reglist[[i]]$onsets)))
    #e1 <- evaluate(reg, gsamples, use_conv=TRUE)
    #e2 <- evaluate(reg, gsamples)
  }
  
  gen_prediction <- function(a, params) {
    Reduce("+", lapply(1:length(reglist), function(i) {
      gen_ith_prediction(i, a[i], params)
    }))
  }
  
  grid <- if (hrf == "gaussian") {
    expand.grid(sd=c(1,2,3), mean=c(4,6,8))
  } else if (hrf == "half_cos") {
    expand.grid(h2=c(4,5,6), h3=c(5,6,7), h4=7)
  }
  
  rsq <- apply(grid, 1, function(p) {
    #browser()
    D <- do.call(cbind, lapply(1:length(reglist), function(i) {
      gen_ith_prediction(i, 1, p)
    }))
    
    #D <- gen_prediction(rep(1,length(reglist)), p)
    lm.0 <- lm(y ~ D)
    b0 <- coef(lm.0)[-1]
    summary(lm.0)$r.squared
  })
  
  best <- which.max(rsq)
  gbest <- grid[best,]
  
  D <- do.call(cbind, lapply(1:length(reglist), function(i) gen_ith_prediction(i,1, gbest)))
  colnames(D) <- names(reglist)
  
  y <- y - mean(y)

  lm.0 <- lm(y ~ D)
  b0 <- coef(lm.0)[-1]
  
  np <- ncol(grid)
  
  calc_error <- function(params) {
    p <- params[1:np]
    amplitude <- params[(np+1):length(params)] # Amplitude parameters for each of the 4 conditions
    
    #browser()
    
    # Generate the predicted fMRI signal by convolving the onset vectors with the HRF model
    # You will need to adapt this part to match your specific data structure and onset vectors
    pred <- Reduce("+", lapply(1:length(reglist), function(i) gen_ith_prediction(i,amplitude[i],p)))
    print(params)
    pred[is.na(pred)] <- 0
    # Calculate the sum of squared differences between the predicted and observed signals
    error <- sum((y - pred)^2)
    cval = cor(y,pred)
    print(error)
    print(cval)
    #return(error)
    return(cval)
  }
  
  # Use optim to minimize the error
  optim_result <- optim(
    par = c(unlist(gbest), b0),
    method="L-BFGS-B",
    fn = calc_error,
    lower=c(bounds$lower, rep(-10, length(reglist))),
    upper=c(bounds$upper, rep(10, length(reglist))),
    control=list(maxit=100, fnscale=-1, parscale=c(rep(1,np),rep(1/sd(y), length(b0))))
  )
  

}

