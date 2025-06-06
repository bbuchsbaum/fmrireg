% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmri_lm_config.R
\name{fmri_lm_control}
\alias{fmri_lm_control}
\title{Configuration for fmri_lm fitting}
\usage{
fmri_lm_control(robust_options = list(), ar_options = list())
}
\arguments{
\item{robust_options}{list of robust fitting options. See Details.}

\item{ar_options}{list of autoregressive modelling options. See Details.}
}
\value{
An object of class \code{fmri_lm_config}.
}
\description{
\code{fmri_lm_control()} creates an \code{fmri_lm_config} object collecting all
options for robust and autoregressive modelling. It validates inputs and
applies defaults so downstream functions receive a single structured list.
}
\details{
\code{robust_options} may contain:
\itemize{
\item \code{type} (\code{FALSE}, "huber", "bisquare")
\item \code{k_huber}
\item \code{c_tukey}
\item \code{max_iter}
\item \code{scale_scope} ("run", "global")
\item \code{reestimate_phi} (logical)
}

\code{ar_options} may contain:
\itemize{
\item \code{struct} ("iid", "ar1", "ar2", "arp")
\item \code{p} (order for "arp")
\item \code{iter_gls} (integer number of GLS iterations)
\item \code{global} (logical, use global phi)
\item \code{voxelwise} (logical)
\item \code{exact_first} (logical)
}
}
