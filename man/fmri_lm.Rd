% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmrilm.R
\name{fmri_lm}
\alias{fmri_lm}
\alias{fmri_lm.formula}
\alias{fmri_lm.fmri_model}
\title{Fit a Linear Regression Model for fMRI Data Analysis}
\usage{
fmri_lm(formula, ...)

\method{fmri_lm}{formula}(
  formula,
  block,
  baseline_model = NULL,
  dataset,
  durations = 0,
  drop_empty = TRUE,
  robust = FALSE,
  robust_options = NULL,
  ar_options = NULL,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10,
  use_fast_path = FALSE,
  progress = FALSE,
  ar_voxelwise = FALSE,
  parallel_voxels = FALSE,
  cor_struct = NULL,
  cor_iter = NULL,
  cor_global = NULL,
  ar1_exact_first = NULL,
  ar_p = NULL,
  robust_psi = NULL,
  robust_max_iter = NULL,
  robust_scale_scope = NULL,
  ...
)

\method{fmri_lm}{fmri_model}(
  fmrimod,
  dataset = NULL,
  robust = FALSE,
  robust_options = NULL,
  ar_options = NULL,
  strategy = c("runwise", "chunkwise"),
  nchunks = 10,
  use_fast_path = FALSE,
  progress = FALSE,
  ar_voxelwise = FALSE,
  parallel_voxels = FALSE,
  cor_struct = NULL,
  cor_iter = NULL,
  cor_global = NULL,
  ar1_exact_first = NULL,
  ar_p = NULL,
  robust_psi = NULL,
  robust_max_iter = NULL,
  robust_scale_scope = NULL,
  ...
)
}
\arguments{
\item{formula}{The model formula for experimental events.}

\item{...}{Additional arguments passed to the chosen method.}

\item{block}{The model formula for block structure.}

\item{baseline_model}{(Optional) A \code{baseline_model} object. Default is \code{NULL}.}

\item{dataset}{An \code{fmri_dataset} object containing the time-series data.}

\item{durations}{A vector of event durations. Default is \code{0}.}

\item{drop_empty}{Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.}

\item{robust}{Logical or character. Either \code{FALSE} (no robust fitting),
\code{TRUE} (use Huber), or one of \code{"huber"} or \code{"bisquare"}. Default is \code{FALSE}.}

\item{robust_options}{List of robust fitting options. See Details.}

\item{ar_options}{List of autoregressive modeling options. See Details.}

\item{strategy}{The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.}

\item{nchunks}{Number of data chunks when strategy is \code{"chunkwise"}. Default is \code{10}.}

\item{use_fast_path}{Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.}

\item{progress}{Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.}
}
\value{
An object of class \code{fmri_lm}.

A fitted linear regression model for fMRI data analysis.
}
\description{
\code{fmri_lm} is a generic for fitting fMRI regression models. The
default interface accepts a model formula and dataset. An
alternative method can be used with a preconstructed
\code{fmri_model} object that already contains the design and data.
}
\details{
\code{robust_options} may contain:
\itemize{
\item \code{type}: Character or logical. Type of robust fitting (\code{FALSE}, \code{"huber"}, \code{"bisquare"})
\item \code{k_huber}: Numeric. Tuning constant for Huber's psi (default: 1.345)
\item \code{c_tukey}: Numeric. Tuning constant for Tukey's bisquare psi (default: 4.685)
\item \code{max_iter}: Integer. Maximum IRLS iterations (default: 2)
\item \code{scale_scope}: Character. Scope for scale estimation (\code{"run"} or \code{"global"})
\item \code{reestimate_phi}: Logical. Whether to re-estimate AR parameters after robust fitting
}

\code{ar_options} may contain:
\itemize{
\item \code{struct}: Character. Correlation structure (\code{"iid"}, \code{"ar1"}, \code{"ar2"}, \code{"arp"})
\item \code{p}: Integer. AR order when \code{struct = "arp"}
\item \code{iter_gls}: Integer. Number of GLS iterations (default: 1)
\item \code{global}: Logical. Use global AR coefficients (default: FALSE)
\item \code{voxelwise}: Logical. Estimate AR parameters voxel-wise (default: FALSE)
\item \code{exact_first}: Logical. Apply exact AR(1) scaling to first sample (default: FALSE)
}
}
\examples{

facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
header=TRUE), face_gen != "n/a")
facedes$face_gen <- droplevels(factor(facedes$face_gen))
sframe <- sampling_frame(rep(430/2,6), TR=2)
ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
block= ~ run, sampling_frame=sframe)
globonsets <- fmrihrf::global_onsets(sframe, facedes$onset, fmrihrf::blockids(ev))
reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=HRF_GAUSSIAN)
reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=HRF_GAUSSIAN)
time <- samples(sframe, global=TRUE)
y1 <- fmrihrf::evaluate(reg1_signal, time)*1.5
y2 <- fmrihrf::evaluate(reg2_signal, time)*3.0
y <- y1+y2
ys1 <- y + rnorm(length(y), sd=.02)
ys2 <- y + rnorm(length(y), sd=.02)

h <<- gen_hrf(hrf_bspline, N=7, span=25)
dset <- matrix_dataset(cbind(ys1,ys2), TR=2, run_length=fmrihrf::blocklens(sframe), event_table=facedes)
flm <- fmri_lm(onset ~ hrf(face_gen, basis=gen_hrf(hrf_bspline, N=7, span=25)), block = ~ run, 
strategy="chunkwise", nchunks=1, dataset=dset)

}
\seealso{
\code{\link{fmri_dataset}}, \code{\link{fmri_lm_fit}}, \code{\link{fmri_lm_control}}
}
