% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/reg-methods.R
\name{durations}
\alias{durations}
\alias{durations.Reg}
\title{Get event durations from an object}
\usage{
durations(x)

\method{durations}{Reg}(x)
}
\arguments{
\item{x}{The object containing event information (typically an event_term or event_model)}
}
\value{
A numeric vector of durations in seconds, where:
\itemize{
\item Each value represents how long an event lasts
\item Zero values indicate instantaneous events
\item Order matches the corresponding event sequence
}
}
\description{
Extract the duration of events from a model object. Durations represent how long each
event lasts in an fMRI design, typically in seconds. These are important for:
\itemize{
\item Modeling block designs where stimuli have non-zero duration
\item Creating accurate HRF convolutions for extended events
\item Distinguishing between brief and sustained neural activity
}
}
\examples{
# Create event data with varying durations
event_data <- data.frame(
  condition = factor(c("block", "event", "block", "event")),
  onsets = c(1, 10, 20, 30),
  durations = c(8, 0, 8, 0),  # 8s blocks and instantaneous events
  run = c(1, 1, 1, 1)
)

# Create event term
eterm <- event_term(
  list(condition = event_data$condition),
  onsets = event_data$onsets,
  durations = event_data$durations,
  blockids = event_data$run
)

# Get durations from term
dur <- durations(eterm)  # Returns: c(8, 0, 8, 0)
}
\seealso{
\code{\link[=onsets]{onsets()}}, \code{\link[=event_term]{event_term()}}

Other timing: 
\code{\link{global_onsets}()},
\code{\link{onsets}()},
\code{\link{samples}()}
}
\concept{timing}
