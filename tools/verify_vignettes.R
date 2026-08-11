#!/usr/bin/env Rscript

# Render every package vignette twice, using a new R process for every render.
# Run from the repository root:
#
#   Rscript tools/verify_vignettes.R /tmp/fmrireg-vignette-gate

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) {
  stop("Run tools/verify_vignettes.R from the fmrireg repository root.")
}
if (!requireNamespace("xml2", quietly = TRUE)) {
  stop("The verification gate requires the 'xml2' package for figure-alt checks.")
}

rendered_warning_count <- function(html_doc) {
  text_lines <- strsplit(xml2::xml_text(html_doc), "\n", fixed = TRUE)[[1L]]
  sum(grepl(
    "^\\s*(?:#>|##)\\s*Warning(?::|\\b)", text_lines, perl = TRUE
  ))
}

# Guard the guard: knitr's collapsed output is HTML-escaped and does not use
# the "## Warning" form that older versions of this verifier searched for.
warning_probe <- xml2::read_html(
  "<html><body><pre><code>#&gt; Warning: verifier probe</code></pre></body></html>"
)
if (!identical(rendered_warning_count(warning_probe), 1L)) {
  stop("Internal verifier error: knitr warning probe was not detected.")
}

output_root <- if (length(args)) args[[1L]] else tempfile("fmrireg-vignette-gate-")
if (dir.exists(output_root) && length(list.files(output_root, all.files = TRUE,
                                                no.. = TRUE))) {
  stop("Output directory already exists and is not empty: ", output_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)

vignettes <- sort(list.files(
  file.path(repo, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE
))
if (length(vignettes) != 11L) {
  stop("Expected 11 vignettes, found ", length(vignettes), ".")
}

# A successful knit is not evidence if a chunk or helper silently discards
# warnings. Keep every warning observable to the fresh-process handler below.
source_lines <- lapply(vignettes, readLines, warn = FALSE)
names(source_lines) <- basename(vignettes)
hidden_warning_lines <- unlist(lapply(source_lines, function(x) {
  grep("warning\\s*=\\s*FALSE|suppressWarnings\\s*\\(", x, perl = TRUE)
}))
if (length(hidden_warning_lines)) {
  offenders <- unlist(Map(function(name, lines) {
    hits <- grep("warning\\s*=\\s*FALSE|suppressWarnings\\s*\\(",
                 lines, perl = TRUE)
    if (!length(hits)) return(character())
    paste0(name, ":", hits)
  }, names(source_lines), source_lines), use.names = FALSE)
  stop(
    "Vignette source suppresses warnings at: ",
    paste(offenders, collapse = ", "),
    ". Qualify expected host warnings in the render harness instead."
  )
}

r_literal <- function(x) encodeString(x, quote = "\"")
records <- vector("list", 2L * length(vignettes))
record_i <- 0L

for (pass in 1:2) {
  pass_dir <- file.path(output_root, paste0("pass-", pass))
  dir.create(pass_dir)

  for (input in vignettes) {
    stem <- tools::file_path_sans_ext(basename(input))
    log_file <- file.path(pass_dir, paste0(stem, ".log"))
    expression <- paste0(
      "options(warn = 1); ",
      "withCallingHandlers({ ",
      "pkgload::load_all(", r_literal(repo),
      ", export_all = FALSE, quiet = TRUE); ",
      "rmarkdown::render(", r_literal(input),
      ", output_dir = ", r_literal(pass_dir),
      ", quiet = TRUE, envir = new.env(parent = globalenv())) ",
      "}, ",
      "warning = function(w) { ",
      "msg <- conditionMessage(w); ",
      "if (grepl(\"^package 'testthat' was built under R version\", msg)) { ",
      "message(\"Qualified host note: \" , msg); ",
      "invokeRestart(\"muffleWarning\") ",
      "} else stop(w) })"
    )

    started <- proc.time()[["elapsed"]]
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", "-e", shQuote(expression)),
      stdout = log_file,
      stderr = log_file,
      env = c("LANG=C", "LC_ALL=C")
    )
    elapsed <- proc.time()[["elapsed"]] - started
    html_file <- file.path(pass_dir, paste0(stem, ".html"))
    log_lines <- readLines(log_file, warn = FALSE)
    html_doc <- if (file.exists(html_file)) {
      xml2::read_html(html_file)
    } else {
      NULL
    }
    figure_alts <- if (!is.null(html_doc)) {
      xml2::xml_attr(xml2::xml_find_all(html_doc, ".//img"), "alt")
    } else {
      character()
    }

    record_i <- record_i + 1L
    records[[record_i]] <- data.frame(
      pass = pass,
      vignette = basename(input),
      seconds = unname(elapsed),
      status = status,
      html_md5 = if (file.exists(html_file)) unname(tools::md5sum(html_file)) else NA_character_,
      console_warning_lines = sum(grepl("Warning", log_lines, fixed = TRUE)),
      rendered_warning_blocks = if (!is.null(html_doc)) {
        rendered_warning_count(html_doc)
      } else {
        0L
      },
      figures = length(figure_alts),
      missing_figure_alt = sum(is.na(figure_alts) | !nzchar(trimws(figure_alts))),
      generic_figure_alt = sum(grepl("^Figure illustrating", figure_alts), na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    if (!identical(status, 0L) || !file.exists(html_file)) {
      stop("Render failed for ", basename(input), " in pass ", pass,
           ". See ", log_file, ".")
    }
    message(sprintf("pass %d: %-32s %.2fs", pass, basename(input), elapsed))
  }
}

receipt <- do.call(rbind, records)
receipt_file <- file.path(output_root, "render-receipt.csv")
write.csv(receipt, receipt_file, row.names = FALSE)

hashes <- split(receipt$html_md5, receipt$vignette)
nondeterministic <- names(Filter(function(x) length(unique(x)) != 1L, hashes))
if (length(nondeterministic)) {
  stop(
    "Fresh-session renders were not byte-identical for: ",
    paste(nondeterministic, collapse = ", "),
    ". See ", receipt_file, "."
  )
}

if (any(receipt$console_warning_lines > 0L) ||
    any(receipt$rendered_warning_blocks > 0L)) {
  stop("One or more renders contain unexplained warnings; inspect the receipt and logs.")
}
if (any(receipt$missing_figure_alt > 0L) ||
    any(receipt$generic_figure_alt > 0L)) {
  stop("One or more rendered figures lack specific alt text; inspect the receipt.")
}

message("All 11 vignettes rendered twice with identical HTML hashes.")
message("Receipt: ", receipt_file)
