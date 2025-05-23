% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/event_model.R, R/regressor.R
\name{nbasis}
\alias{nbasis}
\alias{nbasis.convolved_term}
\alias{nbasis.hrfspec}
\title{Return number of basis functions associated with HRF}
\usage{
nbasis(x)

\method{nbasis}{convolved_term}(x)

\method{nbasis}{hrfspec}(x)
}
\arguments{
\item{x}{The object to query (typically an HRF, hrfspec, or convolved_term)}
}
\value{
An integer indicating the number of basis functions:
\itemize{
\item 1 for canonical HRFs (e.g., SPM gamma)
\item >1 for flexible basis sets (e.g., Fourier, B-spline)
\item For convolved terms: number of basis functions per condition
}
}
\description{
Get the number of basis functions used in a hemodynamic response function (HRF) or
model term. For canonical HRFs (like SPM's canonical HRF), this returns 1. For
flexible basis sets (like Fourier or B-spline bases), this returns the number of
basis functions used to model the response shape.
}
\examples{
# Check basis functions for different HRF types

# Canonical HRF (single basis)
canonical_hrf <- HRF_SPMG1
nbasis(canonical_hrf)  # Returns: 1

# Fourier basis set
fourier_hrf <- getHRF("fourier", nbasis = 3)
nbasis(fourier_hrf)  # Returns: 3

# Create event model with multiple basis functions
event_data <- data.frame(
  condition = factor(c("A", "B", "A", "B")),
  onsets = c(1, 10, 20, 30),
  run = c(1, 1, 1, 1)
)
sframe <- sampling_frame(blocklens = 50, TR = 2)

# Model with Fourier basis
evmodel <- event_model(
  onsets ~ hrf(condition, basis = "fourier", nbasis = 3),
  data = event_data,
  block = ~run,
  sampling_frame = sframe
)

# Get number of basis functions for model term
nbasis(evmodel)  # Returns: 3 (basis functions per condition)
}
\seealso{
\code{\link[=HRF_SPMG1]{HRF_SPMG1()}}, \code{\link[=event_model]{event_model()}}

Other hrf: 
\code{\link{fitted_hrf}()}
}
\concept{hrf}
