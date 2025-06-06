% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{predict.ParametricBasis}
\alias{predict.ParametricBasis}
\alias{predict.Standardized}
\alias{predict.Poly}
\alias{predict.BSpline}
\alias{predict.Ident}
\alias{predict.Scale}
\alias{predict.ScaleWithin}
\alias{predict.RobustScale}
\title{Predict from a ParametricBasis object}
\usage{
predict.ParametricBasis(object, newdata, ...)

\method{predict}{Standardized}(object, newdata, ...)

\method{predict}{Poly}(object, newdata, ...)

\method{predict}{BSpline}(object, newdata, ...)

\method{predict}{Ident}(object, newdata, ...)

\method{predict}{Scale}(object, newdata, ...)

\method{predict}{ScaleWithin}(object, newdata, newgroup, ...)

\method{predict}{RobustScale}(object, newdata, ...)
}
\arguments{
\item{object}{ParametricBasis object}

\item{newdata}{Numeric vector of new values}

\item{...}{Additional arguments}

\item{newgroup}{Optional factor for group-dependent bases.}
}
\description{
Dispatch to the appropriate method for transforming new data
according to a specific parametric basis.
}
