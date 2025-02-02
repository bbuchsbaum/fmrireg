% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/event_model.R
\name{contrast_weights}
\alias{contrast_weights}
\alias{contrast_weights.convolved_term}
\alias{contrast_weights.fmri_model}
\alias{contrast_weights.event_model}
\title{contrast_weights}
\usage{
contrast_weights(x, ...)

\method{contrast_weights}{convolved_term}(x, ...)

\method{contrast_weights}{fmri_model}(x, ...)

\method{contrast_weights}{event_model}(x, ...)
}
\arguments{
\item{x}{the object}

\item{...}{extra args}
}
\description{
This function calculates the contrast weights based on the contrast specification
provided by the user. It is a generic function that dispatches to the appropriate
method depending on the class of the contrast
}
\details{
Calculate contrast weights for a given contrast specification and term.
}
