% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afni.R
\name{gen_afni_lm.fmri_config}
\alias{gen_afni_lm.fmri_config}
\title{Generate an AFNI Linear Model from an fmri_config object}
\usage{
\method{gen_afni_lm}{fmri_config}(x, ...)
}
\arguments{
\item{x}{An fmri_config object containing the necessary configuration for the analysis}

\item{...}{Additional arguments passed to the function}
}
\value{
An afni_lm object representing the fitted linear model
}
\description{
This function takes an fmri_config object and generates an AFNI linear model
by reading the necessary files, creating an fmri_dataset, event_model,
baseline_model, and fmri_model, and fitting the model using afni_lm.
}
\keyword{internal}
