% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ar_whiten_inplace}
\alias{ar_whiten_inplace}
\title{AR(p) whitening of data and design matrices}
\usage{
ar_whiten_inplace(Y, X, phi_coeffs, exact_first_ar1 = FALSE, parallel = TRUE)
}
\arguments{
\item{Y}{Numeric matrix of data (time x voxels)}

\item{X}{Numeric matrix of design (time x predictors)}

\item{phi_coeffs}{Numeric vector of AR coefficients (length p)}

\item{exact_first_ar1}{Logical, apply exact variance-normalizing scaling
of first sample for AR(1). For p > 1, no scaling is applied.}

\item{parallel}{Logical, enable OpenMP parallelization across columns}
}
\value{
List with components 'Y' and 'X' containing the whitened matrices
}
\description{
Applies a causal AR filter defined by \code{phi_coeffs} to both \code{Y} and \code{X}
matrices in place. The filter equation is:
}
\details{
v_t = y_t - sum(phi_k * y_{t-k}, k=1 to p)

The function assumes valid (stationary) AR coefficients are provided.
No checks for stationarity are performed.

For exact_first_ar1 = TRUE and p = 1, the first residual is multiplied
by sqrt(1 - phi^2) for proper variance normalization. This scaling is
only applied for AR(1) models.
}
\keyword{internal}
