#' Run the fmrireg command line interface
#'
#' @param args Character vector of command line arguments.
#'
#' @return Integer exit status. Returns `0` on success, `1` for
#'   command-level validation failures, and `2` for usage or runtime errors.
#' @export
fmrireg_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  tryCatch(
    .fmrireg_cli_dispatch(args),
    fmrireg_cli_usage_error = function(err) {
      .cli_stderr(conditionMessage(err))
      2L
    },
    fmrireg_cli_domain_error = function(err) {
      .cli_stderr(conditionMessage(err))
      1L
    },
    error = function(err) {
      .cli_stderr(conditionMessage(err))
      2L
    }
  )
}

#' Install the fmrireg command wrapper
#'
#' @param dest_dir Destination directory for installed command wrappers.
#' @param overwrite Logical; overwrite an existing installed wrapper.
#' @param commands Character vector of command names to install. Defaults to all
#'   available commands.
#'
#' @return Character vector of installed paths, returned invisibly.
#' @export
install_cli <- function(dest_dir = "~/.local/bin",
                        overwrite = FALSE,
                        commands = NULL) {
  exec_dir <- .fmrireg_exec_dir()
  command_map <- c(fmrireg = file.path(exec_dir, "fmrireg"))

  if (is.null(commands)) {
    commands <- names(command_map)
  }

  unknown <- setdiff(commands, names(command_map))
  if (length(unknown) > 0) {
    stop(
      "Unknown command(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  dest_dir <- path.expand(dest_dir)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }

  installed <- character(length(commands))
  for (i in seq_along(commands)) {
    command <- commands[[i]]
    src <- command_map[[command]]
    dest <- file.path(dest_dir, command)

    if (file.exists(dest) && !isTRUE(overwrite)) {
      stop(
        "Refusing to overwrite existing file: ",
        dest,
        ". Set overwrite = TRUE to replace it.",
        call. = FALSE
      )
    }

    ok <- file.copy(src, dest, overwrite = isTRUE(overwrite))
    if (!ok) {
      stop("Failed to copy CLI wrapper to: ", dest, call. = FALSE)
    }

    if (.Platform$OS.type != "windows") {
      Sys.chmod(dest, mode = "755")
    }

    installed[[i]] <- normalizePath(dest, winslash = "/", mustWork = FALSE)
  }

  path_entries <- strsplit(Sys.getenv("PATH", ""), .Platform$path.sep, fixed = TRUE)[[1]]
  normalized_dest <- normalizePath(dest_dir, winslash = "/", mustWork = FALSE)
  normalized_path <- vapply(
    path_entries[nzchar(path_entries)],
    normalizePath,
    character(1),
    winslash = "/",
    mustWork = FALSE
  )

  if (!normalized_dest %in% normalized_path) {
    message(
      "The directory '", normalized_dest,
      "' is not currently on PATH. Add it to use the installed command directly."
    )
  }

  invisible(installed)
}

.fmrireg_cli_dispatch <- function(args) {
  args <- as.character(args %||% character())

  if (length(args) == 0L) {
    .print_top_level_help()
    return(0L)
  }

  first <- args[[1]]

  if (first %in% c("-h", "--help")) {
    .print_top_level_help()
    return(0L)
  }

  if (identical(first, "help")) {
    return(.handle_help_command(args[-1]))
  }

  if (identical(first, "benchmark")) {
    return(.handle_benchmark_command(args[-1]))
  }

  .usage_error(
    "Unknown command: ", first,
    ". Run 'fmrireg --help' for available commands."
  )
}

.handle_help_command <- function(args) {
  if (length(args) == 0L) {
    .print_top_level_help()
    return(0L)
  }

  topic <- args[[1]]
  if (identical(topic, "benchmark")) {
    .print_benchmark_help()
    return(0L)
  }

  .usage_error("Unknown help topic: ", topic)
}

.handle_benchmark_command <- function(args) {
  if (length(args) == 0L) {
    .print_benchmark_help()
    return(0L)
  }

  subcommand <- args[[1]]
  sub_args <- args[-1]

  if (subcommand %in% c("-h", "--help", "help")) {
    .print_benchmark_help()
    return(0L)
  }

  if (identical(subcommand, "list")) {
    opts <- .parse_cli_options(
      sub_args,
      spec = list(
        help = list(type = "flag", default = FALSE),
        json = list(type = "flag", default = FALSE)
      )
    )
    if (isTRUE(opts$options$help)) {
      .print_benchmark_list_help()
      return(0L)
    }
    return(.run_benchmark_list(json = isTRUE(opts$options$json)))
  }

  if (identical(subcommand, "summary")) {
    opts <- .parse_cli_options(
      sub_args,
      spec = list(
        help = list(type = "flag", default = FALSE),
        json = list(type = "flag", default = FALSE),
        dataset = list(type = "value", required = TRUE)
      )
    )
    if (isTRUE(opts$options$help)) {
      .print_benchmark_summary_help()
      return(0L)
    }
    return(.run_benchmark_summary(
      dataset = opts$options$dataset,
      json = isTRUE(opts$options$json)
    ))
  }

  if (identical(subcommand, "metadata")) {
    opts <- .parse_cli_options(
      sub_args,
      spec = list(
        help = list(type = "flag", default = FALSE),
        json = list(type = "flag", default = FALSE)
      )
    )
    if (isTRUE(opts$options$help)) {
      .print_benchmark_metadata_help()
      return(0L)
    }
    return(.run_benchmark_metadata(json = isTRUE(opts$options$json)))
  }

  .usage_error(
    "Unknown benchmark subcommand: ", subcommand,
    ". Run 'fmrireg benchmark --help' for available subcommands."
  )
}

.run_benchmark_list <- function(json = FALSE) {
  datasets <- tryCatch(
    list_benchmark_datasets(),
    error = function(err) .domain_error(conditionMessage(err))
  )

  if (isTRUE(json)) {
    .write_json(datasets)
  } else {
    .write_table(datasets)
  }

  0L
}

.run_benchmark_summary <- function(dataset, json = FALSE) {
  summary_info <- tryCatch(
    get_benchmark_summary(dataset),
    error = function(err) .domain_error(conditionMessage(err))
  )

  if (isTRUE(json)) {
    .write_json(summary_info)
  } else {
    .write_benchmark_summary(summary_info, dataset)
  }

  0L
}

.run_benchmark_metadata <- function(json = FALSE) {
  metadata <- tryCatch(
    load_benchmark_dataset("metadata"),
    error = function(err) .domain_error(conditionMessage(err))
  )

  if (isTRUE(json)) {
    .write_json(metadata)
  } else {
    .write_named_list(metadata, title = "Benchmark metadata")
  }

  0L
}

.parse_cli_options <- function(args, spec, min_positionals = 0L, max_positionals = 0L) {
  opts <- lapply(spec, function(entry) entry$default %||% NULL)
  positionals <- character()
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--")) {
      positionals <- c(positionals, args[seq.int(i + 1L, length(args))])
      break
    }

    if (!startsWith(arg, "--")) {
      positionals <- c(positionals, arg)
      i <- i + 1L
      next
    }

    if (identical(arg, "--help")) {
      if (!"help" %in% names(spec)) {
        .usage_error("Unsupported option: --help")
      }
      opts$help <- TRUE
      i <- i + 1L
      next
    }

    raw <- substring(arg, 3L)
    negative <- FALSE
    if (startsWith(raw, "no-")) {
      negative <- TRUE
      raw <- substring(raw, 4L)
    }

    pieces <- strsplit(raw, "=", fixed = TRUE)[[1]]
    opt_name <- pieces[[1]]
    opt_value <- NULL
    if (length(pieces) > 1L) {
      opt_value <- paste(pieces[-1], collapse = "=")
    }

    entry <- spec[[opt_name]]
    if (is.null(entry)) {
      .usage_error("Unknown option: --", opt_name)
    }

    if (identical(entry$type, "flag")) {
      if (!is.null(opt_value)) {
        .usage_error("Flag option '--", opt_name, "' does not take a value.")
      }
      opts[[opt_name]] <- !negative
      i <- i + 1L
      next
    }

    if (negative) {
      .usage_error("Option '--", opt_name, "' does not support the --no- form.")
    }

    if (is.null(opt_value)) {
      if (i == length(args)) {
        .usage_error("Missing value for option '--", opt_name, "'.")
      }
      opt_value <- args[[i + 1L]]
      i <- i + 1L
    }

    opts[[opt_name]] <- opt_value
    i <- i + 1L
  }

  if (length(positionals) > 0L) {
    if (length(positionals) < min_positionals || length(positionals) > max_positionals) {
      .usage_error(
        "Unexpected positional argument(s): ",
        paste(positionals, collapse = ", ")
      )
    }
  }

  required <- names(spec)[vapply(spec, function(entry) isTRUE(entry$required), logical(1))]
  missing_required <- required[vapply(required, function(name) {
    is.null(opts[[name]]) || identical(opts[[name]], "")
  }, logical(1))]

  if (length(missing_required) > 0L && !isTRUE(opts$help %||% FALSE)) {
    .usage_error(
      "Missing required option(s): ",
      paste(paste0("--", missing_required), collapse = ", ")
    )
  }

  list(options = opts, positionals = positionals)
}

.write_json <- function(x) {
  cat(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n", sep = "")
}

.write_table <- function(x) {
  utils::write.table(
    x,
    file = stdout(),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

.write_named_list <- function(x, title = NULL) {
  if (!is.null(title)) {
    cat(title, "\n", sep = "")
  }

  for (nm in names(x)) {
    value_text <- paste(capture.output(str(x[[nm]], give.attr = FALSE)), collapse = " ")
    cat(nm, ": ", value_text, "\n", sep = "")
  }
}

.write_benchmark_summary <- function(summary_info, dataset) {
  cat("Dataset: ", dataset, "\n", sep = "")
  cat("Description: ", summary_info$description, "\n", sep = "")
  cat("Timepoints: ", summary_info$dimensions$n_timepoints, "\n", sep = "")
  cat("Voxels: ", summary_info$dimensions$n_voxels, "\n", sep = "")
  cat("Events: ", summary_info$dimensions$n_events, "\n", sep = "")
  cat("Conditions: ", summary_info$dimensions$n_conditions, "\n", sep = "")
  cat("TR: ", summary_info$experimental_design$TR, "\n", sep = "")
  cat("Total time: ", summary_info$experimental_design$total_time, "\n", sep = "")
  cat("Target SNR: ", summary_info$experimental_design$target_snr, "\n", sep = "")
  cat(
    "Condition labels: ",
    paste(summary_info$experimental_design$conditions, collapse = ", "),
    "\n",
    sep = ""
  )
}

.print_top_level_help <- function() {
  cat(
    paste(
      "Usage:",
      "  fmrireg <command> [<args>]",
      "",
      "Commands:",
      "  benchmark   Inspect bundled benchmark datasets",
      "  help        Show help for a command",
      "",
      "Examples:",
      "  fmrireg benchmark list",
      "  fmrireg benchmark summary --dataset BM_Canonical_HighSNR --json",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

.print_benchmark_help <- function() {
  cat(
    paste(
      "Usage:",
      "  fmrireg benchmark <subcommand> [<args>]",
      "",
      "Subcommands:",
      "  list       List benchmark datasets",
      "  summary    Show summary information for one dataset",
      "  metadata   Show package benchmark metadata",
      "",
      "Examples:",
      "  fmrireg benchmark list --json",
      "  fmrireg benchmark summary --dataset BM_Canonical_HighSNR",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

.print_benchmark_list_help <- function() {
  cat(
    paste(
      "Usage:",
      "  fmrireg benchmark list [--json]",
      "",
      "Options:",
      "  --json   Emit machine-readable JSON",
      "  --help   Show this help",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

.print_benchmark_summary_help <- function() {
  cat(
    paste(
      "Usage:",
      "  fmrireg benchmark summary --dataset <name> [--json]",
      "",
      "Options:",
      "  --dataset <name>   Benchmark dataset name",
      "  --json             Emit machine-readable JSON",
      "  --help             Show this help",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

.print_benchmark_metadata_help <- function() {
  cat(
    paste(
      "Usage:",
      "  fmrireg benchmark metadata [--json]",
      "",
      "Options:",
      "  --json   Emit machine-readable JSON",
      "  --help   Show this help",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

.cli_stderr <- function(...) {
  cat(..., "\n", file = stderr(), sep = "")
}

.usage_error <- function(...) {
  err <- simpleError(paste0(...))
  class(err) <- c("fmrireg_cli_usage_error", class(err))
  stop(err)
}

.domain_error <- function(...) {
  err <- simpleError(paste0(...))
  class(err) <- c("fmrireg_cli_domain_error", class(err))
  stop(err)
}

.fmrireg_exec_dir <- function() {
  installed <- system.file("exec", package = "fmrireg")
  if (nzchar(installed)) {
    return(installed)
  }

  namespace_path <- tryCatch(
    getNamespaceInfo(asNamespace("fmrireg"), "path"),
    error = function(err) ""
  )
  candidate <- file.path(namespace_path, "exec")
  if (dir.exists(candidate)) {
    return(candidate)
  }

  stop("Could not locate installed CLI wrappers for fmrireg.", call. = FALSE)
}
