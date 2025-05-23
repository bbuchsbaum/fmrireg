% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event_vector.R
\name{split_onsets.event_term}
\alias{split_onsets.event_term}
\title{Split onsets of an event_term object.}
\usage{
\method{split_onsets}{event_term}(x, sframe, global = FALSE, blocksplit = FALSE, ...)
}
\arguments{
\item{x}{An event_term object.}

\item{sframe}{A data frame representing the sampling frame.}

\item{global}{Logical; if TRUE, use global onsets (default: FALSE).}

\item{blocksplit}{Logical; if TRUE, split onsets by block IDs (default: FALSE).}

\item{...}{Additional arguments.}
}
\value{
A list of numeric vectors for each factor level or block.
}
\description{
Splits the onsets based on factor levels or block IDs.
}
