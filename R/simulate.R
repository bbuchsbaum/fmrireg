
sim_ts <- function(ncond, nreps=12, amps=rep(1,ncond), isi=c(3,6), TR=1.5) {
  cond <- letters[1:ncond]
  trials <- sample(rep(cond, nreps))
  isis <- sample(isi[1]:isi[2], length(trials), replace=TRUE)
  onset <- cumsum(isis)
  
  time <- seq(0, max(onset+12), by=1.5)
  ymat <- do.call(cbind, lapply(1:length(cond), function(i) {
    print(i)
    idx <- which(trials == cond[i])
    reg <- regressor(onset[idx], amplitude=amps[i])
    evaluate(reg, time)
  }))
}