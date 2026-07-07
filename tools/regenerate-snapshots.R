#!/usr/bin/env Rscript

# Regenerates the reference values pinned in the seeded-drift tripwire files
# (inst/tinytest/test-reproducibility-*.R). Each file is parsed, evaluated
# top to bottom against the installed package to recompute its fit(s), and
# every top-level `referenceX <- list(...)` assignment is rewritten with
# freshly computed values, read off the expression each field is checked
# against (`expect_equal(<expr>, referenceX$field)`). Field order is taken
# from the existing list() call, so unrelated formatting is untouched.
#
# Usage: Rscript tools/regenerate-snapshots.R

suppressPackageStartupMessages({
  library(dbarts)
  library(tinytest)
})

scriptPath <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
repoRoot <- normalizePath(file.path(dirname(scriptPath), ".."))

files <- file.path(
  repoRoot,
  "inst",
  "tinytest",
  c(
    "test-reproducibility-continuousResponse-singleThreaded.R",
    "test-reproducibility-continuousResponse-multithreaded.R",
    "test-reproducibility-binaryResponse.R",
    "test-reproducibility-rbart.R",
    "test-reproducibility-xbart.R"
  )
)

isReferenceAssign <- function(e) {
  is.call(e) &&
    identical(e[[1L]], as.name("<-")) &&
    is.name(e[[2L]]) &&
    startsWith(as.character(e[[2L]]), "reference")
}

isFieldAccess <- function(e, refNames) {
  is.call(e) &&
    identical(e[[1L]], as.name("$")) &&
    is.name(e[[2L]]) &&
    as.character(e[[2L]]) %in% refNames
}

regenerateFile <- function(path) {
  lines <- readLines(path)
  exprs <- parse(file = path, keep.source = TRUE)
  srcrefs <- attr(exprs, "srcref")

  refIndex <- Filter(
    function(i) isReferenceAssign(exprs[[i]]),
    seq_along(exprs)
  )
  refNames <- vapply(refIndex, function(i) as.character(exprs[[i]][[2L]]), "")

  # observed-value expression for each reference$field, read off the
  # expect_equal() call that checks it
  fieldExprs <- setNames(vector("list", length(refNames)), refNames)
  for (e in exprs) {
    if (
      is.call(e) &&
        identical(e[[1L]], as.name("expect_equal")) &&
        length(e) >= 3L &&
        isFieldAccess(e[[3L]], refNames)
    ) {
      refName <- as.character(e[[3L]][[2L]])
      fieldName <- as.character(e[[3L]][[3L]])
      fieldExprs[[refName]][[fieldName]] <- e[[2L]]
    }
  }

  # single pass: evaluate expressions in order, but at each reference
  # assignment's original position, substitute freshly computed values
  # instead of re-running it (so later cleanup, e.g. rm(bartFit), doesn't
  # remove what a field expression still needs)
  env <- new.env(parent = globalenv())
  replacements <- list()
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (i %in% refIndex) {
      refName <- as.character(e[[2L]])
      fieldOrder <- names(e[[3L]])[-1L]
      values <- lapply(fieldOrder, function(fn) {
        eval(fieldExprs[[refName]][[fn]], envir = env)
      })
      names(values) <- fieldOrder
      assign(refName, values, envir = env)
      text <- paste0(refName, " <- ", paste(deparse(values), collapse = "\n"))
      replacements[[length(replacements) + 1L]] <- list(
        srcref = srcrefs[[i]],
        text = text
      )
    } else {
      eval(e, envir = env)
    }
  }

  # apply in reverse source order so earlier spans stay valid
  order <- order(
    vapply(replacements, function(r) r$srcref[1L], 0L),
    decreasing = TRUE
  )
  for (r in replacements[order]) {
    sr <- r$srcref
    l1 <- sr[1L]
    c1 <- sr[2L]
    l2 <- sr[3L]
    c2 <- sr[4L]
    prefix <- substr(lines[l1], 1L, c1 - 1L)
    suffix <- substr(lines[l2], c2 + 1L, nchar(lines[l2]))
    spliced <- strsplit(paste0(prefix, r$text, suffix), "\n", fixed = TRUE)[[
      1L
    ]]
    lines <- c(
      if (l1 > 1L) lines[seq_len(l1 - 1L)],
      spliced,
      if (l2 < length(lines)) lines[(l2 + 1L):length(lines)]
    )
  }

  writeLines(lines, path)

  air <- Sys.which("air")
  if (nzchar(air)) {
    system2(air, c("format", path))
  } else {
    warning("air not found on PATH; ", path, " was not reformatted")
  }
}

invisible(lapply(files, regenerateFile))
