# Shadow base::all in globalenv so test-file-level all() calls with a
# zero-length operand are recorded. Package internals resolve all() through
# their own namespace, so only test-file calls are instrumented.
suppressPackageStartupMessages(library(dbarts))
.vac <- new.env()
.vac$hits <- character(0)
all <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (length(args) > 0L && base::all(vapply(args, function(a) length(a) == 0L, NA))) {
    cl <- paste(deparse(sys.call(), width.cutoff = 200L), collapse = " ")
    .vac$hits <- c(.vac$hits, cl)
  }
  base::all(..., na.rm = na.rm)
}
environment(all) <- globalenv()
args <- commandArgs(trailingOnly = TRUE)
files <- readLines(args[1]); root <- Sys.getenv("SRCROOT")
for (f in files) {
  n0 <- length(.vac$hits)
  invisible(tryCatch(
    suppressWarnings(tinytest::run_test_file(file.path(root, f), verbose = 0)),
    error = function(e) NULL))
  if (length(.vac$hits) > n0)
    cat(basename(f), ":", paste(unique(.vac$hits[(n0+1L):length(.vac$hits)]), collapse=" | "), "\n")
}
cat("TOTAL VACUOUS all() CALLS:", length(.vac$hits), "\n")
