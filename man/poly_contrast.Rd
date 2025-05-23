% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast.R
\name{poly_contrast}
\alias{poly_contrast}
\title{Polynomial Contrast}
\usage{
poly_contrast(A, name, where = NULL, degree = 1, value_map = NULL)
}
\arguments{
\item{A}{A formula specifying the ordered factor.}

\item{name}{A character string identifying the contrast.}

\item{where}{An optional formula for subsetting the data.}

\item{degree}{An integer specifying the degree of the polynomial (default: 1).}

\item{value_map}{An optional list mapping factor levels to numeric values.}
}
\value{
A poly_contrast_spec object containing the specification for generating
polynomial contrast weights.
}
\description{
Create polynomial contrasts for testing trends across ordered factor levels. This is
particularly useful for analyzing factors with a natural ordering (e.g., time, dose).
}
\details{
The function creates orthogonal polynomial contrasts up to the specified degree.
These contrasts can test for linear, quadratic, cubic, and higher-order trends
in the data. The value_map parameter allows for non-uniform spacing between levels.
}
\examples{
# Linear trend across time points
pcon <- poly_contrast(~ time, name = "linear_time", degree = 1)

# Cubic trend with custom spacing
pcon <- poly_contrast(~ dose, name = "dose_cubic",
                     degree = 3,
                     value_map = list("low" = 0, "med" = 2, "high" = 5))

}
\seealso{
\code{\link{oneway_contrast}} for categorical contrasts,
\code{\link{interaction_contrast}} for interaction effects
}
