library(purrr)

## generate an experimental design
gen_design <- function(nstim=50, isi_range=c(0,4), noise_sd=.8, rho=.12, TR=2, rnum=1) {
  isi <- runif(nstim, min=isi_range[1]+.1, max=isi_range[2]) 
  onsets <- cumsum(isi) + 1
  runlen <- ceiling(onsets[length(onsets)] + 16)
  nscans <- round(runlen/TR)
  
  fac <- factor(sample(rep(c("a", "b"), each=nstim/2)))
  
  design <- data.frame(
    isi=isi,
    trial=ordered(rep(1:nstim)),
    fac=fac,
    rnum=rnum,
    onsets=onsets
  )
}

## given a design matrix and number oscans, generate an event model
gen_event_model <- function(desmat, nscans, TR=2) {
  sframe <- sampling_frame(nscans, TR)
  desmat$constant <- factor(rep(1, nrow(desmat)))
  emod <- event_model(onsets ~ hrf(constant) + trialwise(fac), 
                      block= ~ rnum, 
                      sampling_frame=sframe,data=desmat)
  
}

## given an experimental design, generate a response
gen_simdata <- function(des, rho=.12, TR=2, noise_sd=.8, rnum=1) {
  runlen <- ceiling(des$onsets[nrow(des)] + 16)
  nscans <- round(runlen/TR)
  
  baseline <- rnorm(nscans, mean=0, sd=noise_sd)
  baseline <- c(baseline[1], map_dbl(2:length(baseline), function(i) {
    rho*baseline[i-1] + (1-rho)*baseline[i]
  }))
  
  amp_a <- rnorm(sum(des$fac == "a"), mean=5, sd=.5)
  y_a <- regressor(des$onsets[des$fac == "a"], amplitude=amp_a)
  amp_b <- rnorm(sum(des$fac == "b"), mean=3, sd=.5)
  y_b <- regressor(des$onsets[des$fac == "b"], amplitude=amp_b)
  
  y_a <- evaluate(y_a, seq(1,runlen, length.out=nscans))
  y_b <- evaluate(y_b, seq(1,runlen, length.out=nscans))
  y_signal <- y_a + y_b
  y <- y_signal + baseline
  data.frame(y_a=y_a,y_b=y_b,y_signal=y_signal,y=y, rnum=rnum)
  
}

## generate a full simulation data set for a single iteration
gen_sim_set <- function(nruns=3, isi_range=c(0,4), noise_sd=.8, TR=2) {
  
  des <- seq(1, nruns) %>% 
    map(~ gen_design(nstim=50, isi_range=isi_range, TR=TR, rnum=.)) 
  
  
  ydat <- des %>% 
    imap(~ gen_simdata(.x,TR=2, noise_sd=noise_sd, rnum=.y)) %>% 
    map_df(dplyr::bind_rows)
  
  desmat <- do.call(rbind, des)
  
  ev <- gen_event_model(desmat, nscans=table(ydat$rnum), TR=2)
  
  list(design=desmat, simdata=ydat, ev=ev)
}

test_that("can generate a trialwise regression model with a summed term term", {
  sim <- gen_sim_set(3, isi_range=c(0,4), noise_sd=.8)
 
})

