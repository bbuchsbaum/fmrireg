% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmrilm.R
\name{fmri_lm}
\alias{fmri_lm}
\title{Fit a Linear Regression Model for fMRI Data Analysis}
\usage{
fmri_lm(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations = 0,
  drop_empty = TRUE,
  robust = FALSE,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10,
  ...
)
}
\arguments{
\item{formula}{The model formula for experimental events.}

\item{block}{The model formula for block structure.}

\item{baseline_model}{(Optional) A \code{baseline_model} object. Default is \code{NULL}.}

\item{dataset}{An \code{fmri_dataset} object containing the time-series data.}

\item{durations}{A vector of event durations. Default is \code{0}.}

\item{drop_empty}{Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.}

\item{robust}{Logical. Whether to use robust fitting. Default is \code{FALSE}.}

\item{strategy}{The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.}

\item{nchunks}{Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.}

\item{...}{Additional arguments.}
}
\value{
A fitted linear regression model for fMRI data analysis.
}
\description{
This function fits a linear regression model for fMRI data analysis using the specified model formula,
block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy,
and robust fitting can be enabled if desired.
}
\examples{

facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
header=TRUE), face_gen != "n/a")
facedes$face_gen <- droplevels(factor(facedes$face_gen))
sframe <- sampling_frame(rep(430/2,6), TR=2)
ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
block= ~ run, sampling_frame=sframe)
globonsets <- global_onsets(sframe, facedes$onset, blockids(ev))
reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=HRF_GAUSSIAN)
reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=HRF_GAUSSIAN)
time <- samples(sframe, global=TRUE)
y1 <- evaluate(reg1_signal, time)*1.5
y2 <- evaluate(reg2_signal, time)*3.0
y <- y1+y2
ys1 <- y + rnorm(length(y), sd=.02)
ys2 <- y + rnorm(length(y), sd=.02)

h <<- gen_hrf(hrf_bspline, N=7, span=25)
dset <- matrix_dataset(cbind(ys1,ys2), TR=2, run_length=sframe$blocklens, event_table=facedes)
flm <- fmri_lm(onset ~ hrf(face_gen, basis=gen_hrf(hrf_bspline, N=7, span=25)), block = ~ run, 
strategy="chunkwise", nchunks=1, dataset=dset)

}
\seealso{
\code{\link{fmri_dataset}}, \code{\link{fmri_lm_fit}}
}
