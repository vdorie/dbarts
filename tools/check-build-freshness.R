#!/usr/bin/env Rscript

# Guards against testing a stale install: `R CMD INSTALL` without
# `--preclean` after a header edit can silently keep pre-edit compiled
# objects, so the installed library and the source tree's headers disagree
# while every test still links and runs. Compares the installed dbarts's
# build time against the newest mtime of every hand-authored file under the
# source tree's R/, src/, inst/; exits 1 and names every file that postdates
# the build. Compiled objects/libraries and configure's generated,
# non-".in" outputs are excluded: `./configure` regenerates them before
# compilation finishes, so on a wholly fresh build they land seconds after
# the install's own Built stamp and would otherwise read as false staleness.
#
# Usage: Rscript tools/check-build-freshness.R <lib> <source-tree>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: check-build-freshness.R <lib> <source-tree>")
}
lib <- args[1L]
tree <- args[2L]

built <- packageDescription("dbarts", lib.loc = lib)$Built
if (is.null(built)) {
  stop("'dbarts' at '", lib, "' has no Built field")
}
ts <- sub(" UTC$", "", trimws(strsplit(built, ";")[[1]][3]))
buildTime <- as.POSIXct(ts, tz = "UTC")

generated <- file.path(
  tree,
  c(
    "src/Makevars",
    "src/config.hpp",
    "src/misc/config.h",
    "src/external/config.h",
    "src/include/misc/types.h"
  )
)
files <- list.files(
  file.path(tree, c("R", "src", "inst")),
  recursive = TRUE,
  full.names = TRUE
)
files <- setdiff(files, generated)
files <- files[!grepl("[.](o|a|so|dylib|dll)$", files)]
newer <- files[file.info(files)$mtime > buildTime]
if (length(newer) > 0L) {
  cat(
    "stale install: newer than the build (",
    format(buildTime),
    "):\n",
    sep = ""
  )
  cat(paste0("  ", newer), sep = "\n")
  quit(status = 1L)
}
quit(status = 0L)
