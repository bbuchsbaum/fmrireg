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