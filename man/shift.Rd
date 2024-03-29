% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{shift}
\alias{shift}
\title{Shift a time series object}
\usage{
shift(x, ...)
}
\arguments{
\item{x}{An object representing a time series or a time-based data structure.}

\item{shift_amount}{A numeric value indicating the amount of time to shift the object by.
Positive values will shift the object to the right, while negative values will shift it to the left.}
}
\value{
An object of the same class as the input, shifted by the specified amount.
}
\description{
This is a generic function to shift time series objects. The goal is to
provide a simple way to apply time shifts to various time series objects,
such as regressors and time series data.
}
\examples{
\dontrun{
# Shift a regressor object
shifted_regressor <- shift(my_regressor, 5)

# Shift a time series object
shifted_time_series <- shift(my_time_series, -2)
}

}
