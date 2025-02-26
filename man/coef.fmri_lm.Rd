% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmrilm.R
\name{coef.fmri_lm}
\alias{coef.fmri_lm}
\title{Extract Model Coefficients from an fmri_lm Object}
\usage{
\method{coef}{fmri_lm}(object, type = c("betas", "contrasts"), recon = FALSE, ...)
}
\arguments{
\item{object}{An \code{fmri_lm} object.}

\item{type}{The type of coefficients to extract: \code{"betas"} or \code{"contrasts"}. Default is \code{"betas"}.}

\item{recon}{Logical. If \code{TRUE}, reconstructs the coefficients into a neuroimaging volume. Default is \code{FALSE}.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A tibble or matrix of coefficients.
}
\description{
This function extracts model coefficients (estimates) from an \code{fmri_lm} object.
}
