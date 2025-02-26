% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmridsl_builder.R
\name{create_hrf_object}
\alias{create_hrf_object}
\title{Create HRF Object from Specification}
\usage{
create_hrf_object(hrf_spec)
}
\arguments{
\item{hrf_spec}{A list containing HRF specifications:
\itemize{
\item type: Type of HRF (e.g., "HRF_SPMG1", "HRF_GAUSSIAN")
\item parameters: Optional list of parameters specific to the HRF type
\item definition: Required for custom HRFs, name of function
}}
}
\value{
An HRF object
}
\description{
Create HRF Object from Specification
}
\examples{
# Basic canonical HRF
create_hrf_object(list(type = "HRF_SPMG1"))

# Gaussian HRF with custom parameters
create_hrf_object(list(
  type = "HRF_GAUSSIAN",
  parameters = list(mean = 5, sd = 1)
))

}
