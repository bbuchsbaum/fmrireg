% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afni.R
\name{afni_lm}
\alias{afni_lm}
\title{afni_lm}
\usage{
afni_lm(fmri_mod, dataset, working_dir = ".", polort = -1, jobs = 1,
  censor = NULL, options = list())
}
\arguments{
\item{censor}{}
}
\description{
afni_lm
}
\examples{
etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
dset <- fmri_dataset(scans=c("s1.nii", "s2.nii"), mask="mask.nii", TR=1, run_length=c(50,50),event_table=etab)

emodel <- event_model(onset ~ hrf(fac), block = ~ run, data=etab, sampling_frame=dset$sampling_frame)
bmodel <- baseline_model("bs", degree=4, sframe=dset$sampling_frame)
fmod <- fmri_model(emodel, bmodel)
alm <- afni_lm(fmod, dset)
}