## Helper to attach fmrigds for tests when available
## Ensures public symbols like plan/explain are present without using :::
if (requireNamespace("fmrigds", quietly = TRUE)) {
  suppressPackageStartupMessages(library(fmrigds))
}

