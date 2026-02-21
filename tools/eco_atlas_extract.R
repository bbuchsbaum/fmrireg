#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(digest)
  library(desc)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

slugify <- function(x) {
  x <- tolower(x %||% "untitled")
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-+|-+$)", "", x)
  if (nchar(x) == 0) "untitled" else x
}

write_jsonl <- function(path, records) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  for (rec in records) writeLines(toJSON(rec, auto_unbox = TRUE, null = "null"), con)
}

read_lines <- function(path) {
  if (!file.exists(path)) return(character())
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

get_git_sha <- function() {
  sha <- Sys.getenv("GITHUB_SHA", "")
  if (nzchar(sha)) return(sha)
  sha <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "")
  if (length(sha) == 1) sha else ""
}

parse_namespace_exports <- function(ns_path = "NAMESPACE") {
  if (!file.exists(ns_path)) return(character())
  lines <- read_lines(ns_path)
  exports <- character()

  exp_lines <- grep("^\\s*export\\s*\\(", lines, value = TRUE)
  for (ln in exp_lines) {
    inside <- sub("^\\s*export\\s*\\(", "", ln)
    inside <- sub("\\)\\s*$", "", inside)
    parts <- strsplit(inside, ",")[[1]]
    parts <- trimws(parts)
    parts <- gsub("^['\"]|['\"]$", "", parts)
    exports <- c(exports, parts[nzchar(parts)])
  }

  unique(exports)
}

extract_roxygen_summary <- function(file_lines, func_start_line) {
  if (func_start_line <= 1 || func_start_line > length(file_lines)) return(NULL)

  i <- func_start_line - 1L
  roxy <- character()
  while (i >= 1) {
    line <- file_lines[[i]]
    if (!grepl("^\\s*#'", line)) break
    roxy <- c(roxy, line)
    i <- i - 1L
  }
  if (length(roxy) == 0) return(NULL)
  roxy <- rev(roxy)
  roxy <- trimws(sub("^\\s*#'\\s?", "", roxy))

  candidates <- roxy[nzchar(roxy) & !grepl("^@", roxy)]
  if (length(candidates) == 0) return(NULL)
  candidates[[1]]
}

pairlist_to_signature <- function(fn_call) {
  if (!is.call(fn_call) || as.character(fn_call[[1]]) != "function") return(NULL)
  args <- fn_call[[2]]
  args_str <- paste(deparse(args), collapse = " ")
  args_str <- sub("^pairlist\\(", "", args_str)
  args_str <- sub("\\)$", "", args_str)
  args_str
}

collect_r_defs <- function(r_files) {
  defs <- list()
  for (rf in r_files) {
    file_lines <- read_lines(rf)
    exprs <- tryCatch(parse(rf, keep.source = TRUE), error = function(e) NULL)
    if (is.null(exprs)) next

    for (ex in exprs) {
      if (!is.call(ex)) next
      op_char <- tryCatch(as.character(ex[[1]]), error = function(e) character(0))
      if (length(op_char) != 1) next
      op <- op_char
      if (!(op %in% c("<-", "="))) next
      if (length(ex) < 3) next
      lhs <- ex[[2]]
      rhs <- ex[[3]]

      if (!is.symbol(lhs)) next
      nm <- as.character(lhs)

      if (!is.call(rhs)) next
      fn_head <- tryCatch(as.character(rhs[[1]])[[1]], error = function(e) "")
      if (!identical(fn_head, "function")) next

      sr <- attr(ex, "srcref")
      start_line <- NA_integer_
      end_line <- NA_integer_
      if (!is.null(sr)) {
        start_line <- sr[[1]]
        end_line <- sr[[3]]
      }

      sig_args <- pairlist_to_signature(rhs)
      sig_raw <- if (!is.null(sig_args)) paste0(nm, "(", sig_args, ")") else paste0(nm, "(...)")
      sig <- gsub("as\\.pairlist\\(alist\\(|\\)\\)$", "", sig_raw)
      sig <- gsub("\\)$", ")", sig)

      summary <- if (!is.na(start_line)) extract_roxygen_summary(file_lines, start_line) else NULL

      src <- if (!is.na(start_line) && !is.na(end_line)) {
        list(path = rf, lines = list(start_line, end_line))
      } else {
        NULL
      }

      defs[[nm]] <- list(
        name = nm,
        signature = sig,
        summary = summary,
        source = src
      )
    }
  }
  defs
}

extract_code_fences <- function(path) {
  lines <- read_lines(path)
  out <- list()

  i <- 1L
  while (i <= length(lines)) {
    line <- lines[[i]]

    is_r_fence_start <- grepl("^```\\s*\\{r[^}]*\\}\\s*$", line) || grepl("^```\\s*r\\s*$", line)
    if (!is_r_fence_start) {
      i <- i + 1L
      next
    }

    start_line <- i + 1L
    j <- start_line
    while (j <= length(lines) && !grepl("^```\\s*$", lines[[j]])) j <- j + 1L
    end_line <- j - 1L

    code <- if (end_line >= start_line) paste(lines[start_line:end_line], collapse = "\n") else ""
    out[[length(out) + 1L]] <- list(
      kind = "fence",
      title = basename(path),
      code = code,
      source = list(path = path, lines = c(start_line, end_line))
    )

    i <- j + 1L
  }

  out
}

extract_eco_markers <- function(path) {
  lines <- read_lines(path)
  out <- list()

  for (i in seq_along(lines)) {
    m <- regexec("^\\s*#\\s*ECO:howto\\s+(.*)\\s*$", lines[[i]])
    reg <- regmatches(lines[[i]], m)[[1]]
    if (length(reg) == 0) next
    q <- reg[[2]]

    start <- i + 1L
    end <- min(length(lines), start + 60L)
    for (j in start:end) {
      if (grepl("^\\s*#\\s*ECO:", lines[[j]])) { end <- j - 1L; break }
    }

    code <- if (end >= start) paste(lines[start:end], collapse = "\n") else ""
    out[[length(out) + 1L]] <- list(
      kind = "eco",
      title = q,
      question_seed = q,
      code = code,
      source = list(path = path, lines = c(start, end))
    )
  }

  out
}

extract_test_that_blocks <- function(path) {
  lines <- read_lines(path)
  out <- list()

  count_char <- function(s, ch) {
    n <- gregexpr(ch, s, fixed = TRUE)[[1]]
    if (identical(n, -1L)) 0L else length(n)
  }

  i <- 1L
  while (i <= length(lines)) {
    line <- lines[[i]]
    if (!grepl("test_that\\s*\\(", line)) { i <- i + 1L; next }

    desc <- NA_character_
    m <- regexec("test_that\\s*\\(\\s*[\"']([^\"']+)[\"']", line)
    reg <- regmatches(line, m)[[1]]
    if (length(reg) >= 2) desc <- reg[[2]] else desc <- "test_that"

    j <- i
    started <- FALSE
    depth <- 0L
    captured <- character()

    while (j <= length(lines)) {
      ln <- lines[[j]]

      if (!started && grepl("\\{", ln)) {
        started <- TRUE
      }

      if (started) {
        captured <- c(captured, ln)
      }

      if (started) {
        depth <- depth + count_char(ln, "{") - count_char(ln, "}")
        if (depth <= 0L && grepl("\\}\\s*\\)\\s*$", ln)) {
          break
        }
      }

      if (length(captured) > 120L) break
      j <- j + 1L
    }

    code <- paste(captured, collapse = "\n")

    out[[length(out) + 1L]] <- list(
      kind = "testthat",
      title = desc,
      question_seed = paste0("How do I ", desc, "?"),
      code = code,
      source = list(path = path, lines = c(i, j))
    )

    i <- j + 1L
  }

  out
}

extract_pkg_calls <- function(code) {
  m <- gregexpr("([A-Za-z][A-Za-z0-9._]+)::([A-Za-z][A-Za-z0-9._]+)", code, perl = TRUE)
  hits <- regmatches(code, m)[[1]]
  unique(hits)
}

# -------------------- main --------------------

out_dir <- Sys.getenv("ECO_ATLAS_OUT", "atlas")
max_snips <- as.integer(Sys.getenv("ECO_ATLAS_MAX_SNIPPETS", "200"))
include_tests <- tolower(Sys.getenv("ECO_ATLAS_INCLUDE_TESTS", "0")) %in% c("1", "true", "yes")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- desc::desc(file = "DESCRIPTION")
pkg <- d$get("Package") %||% "unknownpkg"
ver <- d$get("Version") %||% "0.0.0"
sha <- get_git_sha()

manifest <- list(
  package = pkg,
  version = ver,
  language = "R",
  commit_sha = sha,
  build_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  # Keep legacy fields during transition.
  commit = sha,
  built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

exports <- parse_namespace_exports("NAMESPACE")
r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
defs <- collect_r_defs(r_files)

# symbols.jsonl
symbols <- list()
for (nm in exports) {
  def <- defs[[nm]]
  rec <- list(
    symbol = paste0(pkg, "::", nm),
    language = "R",
    type = "function",
    signature = (def$signature %||% paste0(nm, "(...)")),
    summary = def$summary %||% NULL,
    source = def$source %||% list(path = NULL, lines = NULL)
  )
  symbols[[length(symbols) + 1L]] <- rec
}
write_jsonl(file.path(out_dir, "symbols.jsonl"), symbols)

# snippets.jsonl
snippets <- list()

for (rf in r_files) {
  chunks <- extract_eco_markers(rf)
  if (length(chunks) > 0) snippets <- c(snippets, chunks)
}

v_files <- c(
  list.files("vignettes", pattern = "\\.(Rmd|qmd)$", full.names = TRUE),
  list.files("inst", pattern = "\\.(Rmd|qmd)$", full.names = TRUE, recursive = TRUE),
  list.files(".", pattern = "^README\\.(md|Rmd)$", full.names = TRUE)
)
v_files <- unique(v_files[file.exists(v_files)])

for (vf in v_files) {
  chunks <- extract_code_fences(vf)
  if (length(chunks) > 0) {
    for (k in seq_along(chunks)) {
      chunks[[k]]$kind <- if (grepl("^README", basename(vf))) "readme" else "vignette"
      chunks[[k]]$title <- paste0(chunks[[k]]$kind, ": ", chunks[[k]]$title)
    }
    snippets <- c(snippets, chunks)
  }
}

if (include_tests) {
  t_files <- list.files("tests/testthat", pattern = "\\.R$", full.names = TRUE)
  for (tf in t_files) {
    chunks <- extract_test_that_blocks(tf)
    if (length(chunks) > 0) snippets <- c(snippets, chunks)
  }
}

kind_rank <- function(kind) {
  switch(kind, "eco" = 1L, "vignette" = 2L, "readme" = 3L, "testthat" = 4L, 9L)
}

if (length(snippets) > 0) {
  snippets <- snippets[order(vapply(snippets, function(s) kind_rank(s$kind), integer(1)))]
  if (length(snippets) > max_snips) snippets <- snippets[seq_len(max_snips)]
}

snip_records <- list()
edge_records <- list()

for (idx in seq_along(snippets)) {
  s <- snippets[[idx]]
  code <- s$code %||% ""
  title <- s$title %||% "snippet"
  kind <- s$kind %||% "snippet"
  qseed <- s$question_seed %||% NULL

  start_line <- s$source$lines[[1]] %||% NA_integer_
  base_id <- paste0(pkg, "::", kind, "/", slugify(paste0(title, "-", idx, "-", start_line)))

  sym_calls <- extract_pkg_calls(code)
  h <- digest::digest(paste(pkg, kind, title, qseed %||% "", code), algo = "sha1")

  rec <- list(
    id = base_id,
    package = pkg,
    language = "R",
    kind = kind,
    title = title,
    question_seed = qseed,
    code = code,
    symbols = sym_calls,
    source = s$source,
    hash = h
  )

  if (length(sym_calls) > 0) {
    from_symbol <- paste0(pkg, "::snippet_", slugify(paste0(kind, "-", idx, "-", start_line)))

    for (sc in sym_calls) {
      other_pkg <- strsplit(sc, "::", fixed = TRUE)[[1]][1]
      if (!identical(other_pkg, pkg)) {
        edge_records[[length(edge_records) + 1L]] <- list(
          from = from_symbol,
          to = sc,
          source = s$source,
          kind = "uses"
        )
      }
    }
  }

  snip_records[[length(snip_records) + 1L]] <- rec
}

write_jsonl(file.path(out_dir, "snippets.jsonl"), snip_records)
write_jsonl(file.path(out_dir, "edges.jsonl"), edge_records)
cat(sprintf("EcoAtlas extract complete: %s symbols, %s snippets\n", length(symbols), length(snip_records)))
