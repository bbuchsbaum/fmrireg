% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/baseline_model.R,
%   R/event_model.R, R/event_vector.R, R/fmri_model.R
\name{design_matrix}
\alias{design_matrix}
\alias{design_matrix.convolved_term}
\alias{design_matrix.afni_hrf_convolved_term}
\alias{design_matrix.baseline_model}
\alias{design_matrix.baseline_term}
\alias{design_matrix.event_model}
\alias{design_matrix.event_term}
\alias{design_matrix.fmri_model}
\title{design_matrix}
\usage{
design_matrix(x, ...)

\method{design_matrix}{convolved_term}(x, blockid = NULL, ...)

\method{design_matrix}{afni_hrf_convolved_term}(x, blockid = NULL, ...)

\method{design_matrix}{baseline_model}(x, blockid = NULL, allrows = FALSE, ...)

\method{design_matrix}{baseline_term}(x, blockid = NULL, allrows = FALSE, ...)

\method{design_matrix}{event_model}(x, blockid = NULL, ...)

\method{design_matrix}{event_term}(x, drop.empty = TRUE, ...)

\method{design_matrix}{fmri_model}(x, blockid = NULL, ...)
}
\arguments{
\item{x}{The term or model object (typically an event_term, event_model, baseline_model, or fmri_model)}

\item{...}{Additional arguments passed to methods}

\item{blockid}{Numeric vector specifying which blocks/runs to include (optional)}

\item{allrows}{Logical; if TRUE, return all rows including those with no events (default: FALSE)}

\item{drop.empty}{Logical; if TRUE, drop empty conditions (default: TRUE)}
}
\value{
A tibble containing the design matrix, where:
\itemize{
\item Rows represent time points (scans)
\item Columns represent predictor variables
\item Column names indicate the condition or regressor (see Details)
}
}
\description{
Extract or construct the design matrix from a model term or object. The design matrix
contains the predictor variables used in the model, with one row per time point and
one column per predictor. For event-related designs, the design matrix typically
contains the convolved HRF responses. For baseline terms, it contains drift and
nuisance regressors.
}
\details{
\subsection{Column Naming (Post-Refactor)}{

As of version X.Y.Z, column names in the final design matrix generated by
\code{event_model} follow the structure:
\code{term_tag} + \verb{_} + \code{condition_tag} + \verb{_b##} basis suffix

Refer to \code{event_model} documentation for details on \code{term_tag} and \code{condition_tag} generation.
The \code{design_matrix} methods for individual term types (like \code{design_matrix.event_term})
return unconvolved predictors whose names may \emph{not} yet follow this final structure.
The final naming is applied by \code{convolve.event_term} and assembled by
\code{build_event_model_design_matrix}.
}
}
\section{Methods (by class)}{
\itemize{
\item \code{design_matrix(afni_hrf_convolved_term)}: Design matrix for AFNI-convolved terms

AFNI-convolved terms rely on AFNI software for design matrix
construction. This method stops with an informative message.

}}
\seealso{
\code{\link[=event_model]{event_model()}}, \code{\link[=baseline_model]{baseline_model()}}, \code{\link[=fmri_model]{fmri_model()}}
}
\concept{design_matrices}
