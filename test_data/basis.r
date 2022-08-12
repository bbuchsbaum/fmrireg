
out <- do.call(cbind, lapply(1:50000, function(i) {
  nevents <- sample(10:25,1)
  durs = sample(0:10, nevents, replace=TRUE)
  events <- sort(sample(1:100, nevents))
  reg <- regressor(events, duration=durs)
  evaluate(reg, 1:100)
}))

res <- multivarious::pca(out)