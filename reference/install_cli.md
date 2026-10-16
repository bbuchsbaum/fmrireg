# Install the fmrireg command wrapper

Copies the package CLI wrapper from the installed `exec/` directory into
a user-selected destination directory and, on Unix-like systems, marks
the copied file as executable.

## Usage

``` r
install_cli(dest_dir, overwrite = FALSE, commands = NULL)
```

## Arguments

- dest_dir:

  Destination directory for installed command wrappers. Required: there
  is no default, so nothing is written to the home directory unless
  asked (e.g. `install_cli("~/.local/bin")`).

- overwrite:

  Logical; overwrite an existing installed wrapper.

- commands:

  Character vector of command names to install. Defaults to all
  available commands.

## Value

Character vector of installed paths, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
install_cli("~/.local/bin", overwrite = TRUE)
} # }
```
