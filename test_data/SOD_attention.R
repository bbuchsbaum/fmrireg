library(fmrireg)
library(dplyr)
library(mixmeta)
library(multivarious)

atlas <- neuroatlas::get_schaefer_atlas("400", "17")
sids <- (1001:1028)

#hrfb <- gen_hrf(HRF_BSPLINE, N=9)
hrfb <- readRDS("empirical_hrf.rds")
cset <- NULL

elab_freq <- lapply(sids, function(sid) {
  dat <- readRDS(paste0("~/Dropbox/analysis/SOD/fmri/schaefer/", sid, "_schaefer_mean_404.rds"))
  des <- dat$design %>% mutate(elab_onset = elab_onset_ms/1000)
  des$elab <- des$elaboration
  if (sid == "1015" || sid == "1014") {
    des$elab = case_when(
      des$elab == 1 ~ 4,
      des$elab == 2 ~ 3,
      des$elab == 3 ~ 2,
      des$elab == 4 ~ 1, 
      des$elab == 0 ~ 0)
  }
  
  
  des$felab <- factor(des$elab)
  data.frame(subject=sid, elabfreq=table(des$elab), felab=names(table(des$elab)))
  
}) %>% bind_rows()


con1 <- fmrireg::pair_contrast(~ attention == "FA", ~ attention == "DA", name="FA_min_DA_b1", where = ~ basis == "basis01")
con2 <- fmrireg::pair_contrast(~ attention == "FA", ~ attention == "DA", name="FA_min_DA_b2", where = ~ basis == "basis02")
cset <- contrast_set(con1, con2)

run_subject <- function(sid) {
  print(sid)
  dat <- readRDS(paste0("~/Dropbox/analysis/SOD/fmri/schaefer/", sid, "_schaefer_mean_404.rds"))
  des <- dat$design %>% mutate(elab_onset = elab_onset_ms/1000)
  des$elab <- des$elaboration
  if (sid == "1015" || sid == "1014") {
    des$elab = case_when(
                         des$elab == 1 ~ 4,
                         des$elab == 2 ~ 3,
                         des$elab == 3 ~ 2,
                         des$elab == 4 ~ 1, 
                         des$elab == 0 ~ 0)
  }
  
  
  des$felab <- factor(des$elab)
  print(table(des$elab))
  des$elab[des$elab==0] <- NA
  des$elab[is.na(des$elab)] <- mean(des$elab, na.rm=TRUE)
  des$elab <- scale(des$elab)[,1]
  des$accuracy <- factor(ifelse(des$accuracy == 1, "correct", "incorrect"))
  
  dset <- fmrireg::matrix_dataset(as.matrix(dat$data), TR=1.3, run_length=rep(528,4), event_table=des)

  
  ret <- fmrireg::fmri_lm(elab_onset ~ hrf(attention, basis=hrfb, subset=phase=="encoding" & felab!=0, contrasts=cset), #+ 
                            #hrf(attention,elab, basis=hrfb, subset=phase=="encoding"), 
                          strategy="chunkwise",
                        dataset=dset, block = ~ Scan)
  ret
}

res <- lapply(sids, run_subject)
betas <- lapply(res, coef, "contrasts")

out <- lapply(1:404, function(i) {
  b1 <- sapply(1:length(betas), function(j) betas[[j]][i,1])
  b2 <- sapply(1:length(betas), function(j) betas[[j]][i,2])
  
  r1 <- broom::tidy(t.test(unlist(b1))) %>% mutate(basis=1, roi=i)
  r2 <- broom::tidy(t.test(unlist(b2))) %>% mutate(basis=2, roi=i)
  rbind(r1,r2)
  
}) %>% bind_rows()



write.table(out, "~/Dropbox/analysis/SOD/fmri/glm_SOD_attention_only.txt", row.names=FALSE)



