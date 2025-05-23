% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event_model.R
\name{matrix_term}
\alias{matrix_term}
\title{matrix_term}
\usage{
matrix_term(varname, mat)
}
\arguments{
\item{varname}{The name of the variable.}

\item{mat}{The matrix of values.}
}
\value{
A matrix_term object.
}
\description{
Creates a matrix_term object, which is a set of regression variables stored as
a numeric matrix.
}
\details{
Create a matrix_term object.
}
\examples{
mat <- matrix(rnorm(100 * 10), 100, 10)
mterm <- matrix_term("mterm", mat)
}
