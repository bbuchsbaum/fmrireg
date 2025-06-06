% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmrilm.R
\name{reshape_coef}
\alias{reshape_coef}
\title{Reshape Coefficient Data}
\usage{
reshape_coef(df, des, measure = "value")
}
\arguments{
\item{df}{A data frame containing coefficient estimates.}

\item{des}{A data frame containing design information.}

\item{measure}{The name of the value column in the reshaped data. Default is \code{"value"}.}
}
\value{
A data frame in long format with merged design information.
}
\description{
This function reshapes coefficient data from wide to long format and merges it with design information.
}
\keyword{internal}
