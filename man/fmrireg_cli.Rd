% Manual page for CLI entrypoint
\name{fmrireg_cli}
\alias{fmrireg_cli}
\title{Run the fmrireg command line interface}
\usage{
fmrireg_cli(args = commandArgs(trailingOnly = TRUE))
}
\arguments{
\item{args}{Character vector of command line arguments.}
}
\value{
Integer exit status. Returns \code{0} on success, \code{1} for
command-level validation failures, and \code{2} for usage or runtime errors.
}
\description{
Runs the package-owned command line interface for \code{fmrireg}. The current
public command surface focuses on bundled benchmark dataset inspection through
the \code{benchmark} subcommands.
}
\examples{
fmrireg_cli(c("benchmark", "list"))
fmrireg_cli(c(
  "benchmark", "summary",
  "--dataset", "BM_Canonical_HighSNR",
  "--json"
))
}
