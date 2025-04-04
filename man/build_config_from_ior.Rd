% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmridsl.R
\name{build_config_from_ior}
\alias{build_config_from_ior}
\title{Build fmri_config Object from Validated IOR List}
\usage{
build_config_from_ior(validated_ior)
}
\arguments{
\item{validated_ior}{A list returned by \code{parse_and_validate_config}.}
}
\value{
An \code{fmri_config} object ready for model building.
}
\description{
Takes a validated and normalized intermediate representation (IOR) list,
performs context-dependent validation (BIDS checks, event/confound column checks),
infers variable roles, and constructs the final \code{fmri_config} object.
}
\keyword{internal}
