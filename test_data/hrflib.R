onsets = sort(runif(20)*240)

y1 <- evaluate(regressor(onsets, hrf=HRF_GAUSSIAN), seq(1,240,by=1)) + rnorm(240,sd=.01)
y2 <- evaluate(regressor(onsets, hrf=HRF_GAUSSIAN), seq(1,240,by=1)) + rnorm(240,sd=.01)

Y <- cbind(y1,y2)

X <- do.call(cbind, lapply(seq(3,8, by=.1), function(m) {
  evaluate(regressor(onsets, hrf=gen_hrf(HRF_GAUSSIAN, mean=m)), seq(1,240,by=1))
}))

librrary(glmnet)

gres <- glmnet(X,Y, family="mgaussian", pmax=1, nlambda=200, intercept=FALSE)