% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmrilm.R
\name{fmri_lm}
\alias{fmri_lm}
\title{fmri_lm}
\usage{
fmri_lm(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations,
  drop_empty = TRUE,
  robust = FALSE,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10
)
}
\arguments{
\item{formula}{the model formula for experimental events}

\item{block}{the model formula for block structure}

\item{baseline_model}{the \code{baseline_model} object}

\item{dataset}{an object derived from \code{fmri_dataset} containing the time-series data}

\item{durations}{a vector of event durations}

\item{drop_empty}{whether to remove factor levels with size of zero}

\item{robust}{whether to use robust fitting (TRUE or FALSE)}

\item{strategy}{the data splitting strategy}

\item{contrasts}{a set of contrasts}
}
\description{
fmri_lm
}
\examples{
etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
etab2 <- data.frame(onset=c(1,30,65,75), fac=factor(c("A", "B", "A", "B")), run=c(1,1,1,1))
mat <- matrix(rnorm(100*100), 100,100)
dset <- matrix_dataset(mat, TR=1, run_length=c(50,50),event_table=etab)
dset2 <- matrix_dataset(mat, TR=1, run_length=c(100),event_table=etab2)
con  <- pair_contrast(~ fac == "A", ~ fac == "B", name="A_min_B")
con2  <- unit_contrast(~ fac, name="A_min_baseline")
lm.1 <- fmri_lm(onset ~ hrf(fac, contrasts=contrast_set(con,con2)), block= ~ run, dataset=dset)
lm.2 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2)
lm.2a <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, robust=TRUE, dataset=dset2)
lm.3 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2, strategy="runwise")
lm.3a <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, robust=TRUE, dataset=dset2, strategy="runwise")
}
