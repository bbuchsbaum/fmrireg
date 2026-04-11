% Manual page for CLI installer
\name{install_cli}
\alias{install_cli}
\title{Install the fmrireg command wrapper}
\usage{
install_cli(dest_dir = "~/.local/bin", overwrite = FALSE, commands = NULL)
}
\arguments{
\item{dest_dir}{Destination directory for installed command wrappers.}

\item{overwrite}{Logical; overwrite an existing installed wrapper.}

\item{commands}{Character vector of command names to install. Defaults to all
available commands.}
}
\value{
Character vector of installed paths, returned invisibly.
}
\description{
Copies the package CLI wrapper from the installed \code{exec/} directory into a
user-selected destination directory and, on Unix-like systems, marks the copied
file as executable.
}
\examples{
\dontrun{
install_cli("~/.local/bin", overwrite = TRUE)
}
}
