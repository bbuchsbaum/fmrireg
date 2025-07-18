% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{amplitudes}
\alias{amplitudes}
\title{Get event amplitudes from an object}
\usage{
amplitudes(x)
}
\arguments{
\item{x}{The object containing event information (typically an event_term or event_model)}
}
\value{
A numeric vector of amplitude values, where:
\itemize{
\item Each value represents the intensity of an event
\item Default value of 1 indicates unmodulated events
\item Order matches the corresponding event sequence
}
}
\description{
Extract the amplitude or intensity values associated with each event. Amplitudes
represent the strength or magnitude of events and can be used to:
\itemize{
\item Model parametric modulation of neural responses
\item Weight events by their intensity or importance
\item Create amplitude-modulated regressors
}
}
\examples{
# Create event data with varying amplitudes
event_data <- data.frame(
  condition = factor(c("stim", "stim", "stim", "stim")),
  onsets = c(1, 10, 20, 30),
  intensity = c(0.5, 1.0, 1.5, 2.0),  # Parametrically varying intensity
  run = c(1, 1, 1, 1)
)

# Create event term with amplitudes
eterm <- event_term(
  list(condition = event_data$condition),
  onsets = event_data$onsets,
  amplitudes = event_data$intensity,
  blockids = event_data$run
)

# Get amplitudes from term
amp <- amplitudes(eterm)  # Returns: c(0.5, 1.0, 1.5, 2.0)
}
\seealso{
\code{\link[=event_term]{event_term()}}, \code{\link[=onsets]{onsets()}}
}
\concept{event_properties}
