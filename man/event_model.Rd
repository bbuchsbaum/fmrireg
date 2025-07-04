% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/event_model.R
\name{event_model}
\alias{event_model}
\title{Construct an event model}
\usage{
event_model(
  formula_or_list,
  data,
  block,
  sampling_frame,
  durations = 0,
  drop_empty = TRUE,
  precision = 0.3,
  parallel = FALSE,
  progress = FALSE,
  ...
)

event_model(
  formula_or_list,
  data,
  block,
  sampling_frame,
  durations = 0,
  drop_empty = TRUE,
  precision = 0.3,
  parallel = FALSE,
  progress = FALSE,
  ...
)
}
\arguments{
\item{formula_or_list}{Either a formula (e.g., \code{onset ~ hrf(cond) + hrf(mod)})
or a list of pre-defined \code{hrfspec} objects.}

\item{data}{A \code{data.frame} containing event variables referenced in the formula
or needed by the \code{hrfspec} objects.}

\item{block}{A formula (e.g., \code{~ run}) or vector specifying the block/run for each event.}

\item{sampling_frame}{An object of class \code{sampling_frame} defining the scan timing (TR, block lengths).}

\item{durations}{Numeric vector or scalar specifying event durations (seconds). Default is 0.}

\item{drop_empty}{Logical indicating whether to drop empty events during term construction. Default is TRUE.}

\item{precision}{Numeric precision for HRF sampling/convolution. Default is 0.3.}

\item{parallel}{Logical indicating whether to use parallel processing for term convolution (requires \code{future.apply}). Default is FALSE.}

\item{progress}{Logical indicating whether to show a progress bar during term realisation. Default is FALSE.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A list containing the following elements:
\describe{
\item{formula}{The formula used to create the model.}
\item{design}{The design matrix for the model, with one row per time point and one column per predictor variable.}
\item{block_indices}{A list of indices defining the start and end time points of each block.}
\item{baseline}{A vector containing the estimated baseline fMRI signal level for each block.}
\item{dur}{A vector containing the duration (in seconds) of each event or block in the design.}
}

An object of class \code{c("event_model", "list")} containing the terms, design matrix,
sampling frame, and other metadata.
}
\description{
This function creates an event-based fMRI regression model, represented as a data structure.

This is the main constructor for \code{event_model} objects. It unifies the previous
formula and list-based interfaces and uses a more efficient internal pipeline.
}
\details{
\subsection{Column Naming}{

The columns in the resulting design matrix follow the naming convention:
\code{term_tag} + \verb{_} + \code{condition_tag} + \verb{_b##} basis suffix

Where:
\itemize{
\item \code{term_tag}: The unique tag assigned to the \code{hrf()} term (see below).
\item \code{condition_tag}: Represents the specific factor level or continuous regressor within the term (e.g., \code{condition.A}, \code{poly_RT_01}, \code{condition.A_task.go}).
\item \verb{_b##}: Optional suffix added for HRFs with multiple basis functions (e.g., \verb{_b01}, \verb{_b02}).
}
}

\subsection{Term Naming and Clash Resolution}{

Each term in the model (typically defined by an \code{hrf()} call in a formula) gets a
unique \code{term_tag}. This tag is used as the prefix for all columns generated by that term.
\itemize{
\item \strong{Default Naming:} If no explicit \code{id} (or \code{name}) is given in \code{hrf()}, the tag is derived
from the variable names (e.g., \code{hrf(condition)} -> \code{condition}, \code{hrf(RT, acc)} -> \code{RT_acc}).
\item \strong{Explicit Naming:} Use \verb{id=} within \code{hrf()} for an explicit tag (e.g., \code{hrf(condition, id="CondMain")}).
\item \strong{Sanitization:} Dots (\code{.}) in tags are converted to underscores (\verb{_}).
\item \strong{Clash Resolution:} If multiple terms generate the same tag, \verb{#} and a number are appended
to ensure uniqueness (e.g., \code{condition}, \code{condition#1}).
}

This consistent naming scheme replaces the previous \code{compact} and \code{qualified} styles.
}
}
\examples{
# Create a data frame with experimental design
event_data <- data.frame(fac=c("a", "B", "A", "B"), onsets=c(1,10,20,80), run=c(1,1,1,1))

# Create a sampling frame with 50-second blocks and a TR of 2 seconds
sframe <- sampling_frame(blocklens=50, TR=2)

# Create an event model using the `onsets` variable as a predictor, 
#  with a separate baseline for each run
evmodel <- event_model(onsets ~ hrf(onsets), data=event_data, block=~run, sampling_frame=sframe)
dmat <- design_matrix(evmodel)
# Example using formula interface
des <- data.frame(onset = seq(0, 90, by=10),
                  run = rep(1:2, each=5),
                  cond = factor(rep(c("A","B"), 5)),
                  mod = rnorm(10))
sframe <- fmrihrf::sampling_frame(blocklens=c(50, 60), TR=2)

ev_model_form <- event_model(onset ~ hrf(cond) + hrf(mod, basis="spmg3"), 
                            data = des, block = ~run, sampling_frame = sframe)
print(ev_model_form)
head(design_matrix(ev_model_form))

# Example using list interface (less common)
# spec1 <- hrf(cond)
# spec2 <- hrf(mod, basis="spmg3")
# ev_model_list <- event_model(list(spec1, spec2), data=des, block=des$run, sampling_frame=sframe)
# print(ev_model_list)                         
                             
}
