% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naming-utils.R
\name{basis_suffix}
\alias{basis_suffix}
\title{Create Basis Function Suffix}
\usage{
basis_suffix(j, nb)
}
\arguments{
\item{j}{Integer vector of basis indices (1-based).}

\item{nb}{Total number of basis functions.}
}
\value{
Character vector of suffixes (e.g., \verb{_b01}, \verb{_b02}).
}
\description{
Generates the \verb{_b##} suffix for HRF basis functions.
}
\examples{
basis_suffix(1:3, 5)
basis_suffix(1:10, 10)
}
