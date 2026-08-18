#!/usr/bin/env Rscript

# Guards the living design/plan references against silent drift.
#
# 1. INDEX COMPLETENESS: every docs/design/*.md file is listed in
#    docs/design/INDEX.md, and every name INDEX.md lists exists on disk
#    (no phantom rows). docs/plans/INDEX.md and docs/plans/*.md get the
#    same check when the plans index exists.
#
# 2. FEATURE-MATRIX ANCHOR STALENESS: every file:line (and the handful of
#    file(symbol)) anchor docs/design/feature-matrix.md cites, using the
#    ALIAS table at that file's own head, resolves against the live tree -
#    the target file exists, the cited line is within the file's current
#    length, and a symbol quoted next to the anchor is still present
#    somewhere in the target file. Also recomputes the equivalence
#    baseline scenario counts a footnote states (e.g. [f39]) against the
#    actual scenario names in the harness each baseline cites
#    (benchmarks/R/equivalence.R, bcf-equivalence.R,
#    multinomial-equivalence.R) and flags a mismatch.
#
# This is a DEAD-REFERENCE and DRIFTED-COUNT detector, not a semantic
# reviewer: it never re-adjudicates a SHIPPED/REFUSED/etc. cell VALUE, and
# a symbol check only asks whether the string still occurs in the file
# ANYWHERE, not at the cited line - a line number that has merely moved
# within a file that still exists and still contains the symbol is not
# flagged (that is the anchor-resync pass's job). An anchor this guard's
# alias/basename resolution cannot place is skipped, not failed.
#
# Usage: Rscript tools/check-doc-freshness.R [pkg-root]
# Exit status: 0 if every check passes, 1 if any fails (failures print
# with a leading "FAIL:" so they are easy to grep).

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[1L] else "."
p <- function(...) file.path(root, ...)

find <- character(0L)
report <- function(msg) find <<- c(find, msg)

# ---------------------------------------------------------------------
# Part 1: INDEX completeness
# ---------------------------------------------------------------------

# First-column ".md" entries of a "| file | ... |" table, wherever such a
# row appears in the index (multiple tables/sections are all scanned).
indexedNames <- function(indexFile) {
  lines <- readLines(indexFile, warn = FALSE)
  hit <- grep("^\\| *[A-Za-z0-9_.+-]+\\.md *\\|", lines, value = TRUE)
  sub("^\\| *([A-Za-z0-9_.+-]+\\.md) *\\|.*$", "\\1", hit)
}

checkIndex <- function(dir, indexRel) {
  indexFile <- p(indexRel)
  onDisk <- list.files(p(dir), pattern = "\\.md$")
  onDisk <- setdiff(onDisk, basename(indexRel))
  indexed <- unique(indexedNames(indexFile))
  missing <- sort(setdiff(onDisk, indexed))
  phantom <- sort(setdiff(indexed, onDisk))
  for (f in missing) {
    report(sprintf(
      "INDEX: %s/%s exists on disk but is not listed in %s",
      dir,
      f,
      indexRel
    ))
  }
  for (f in phantom) {
    report(sprintf(
      "INDEX: %s/%s is listed in %s but does not exist on disk",
      dir,
      f,
      indexRel
    ))
  }
  length(onDisk)
}

nDesign <- 0L
nPlans <- 0L
if (file.exists(p("docs/design/INDEX.md"))) {
  nDesign <- checkIndex("docs/design", "docs/design/INDEX.md")
} else {
  report("INDEX: docs/design/INDEX.md is missing")
}
if (file.exists(p("docs/plans/INDEX.md"))) {
  nPlans <- checkIndex("docs/plans", "docs/plans/INDEX.md")
}

# ---------------------------------------------------------------------
# Part 2: feature-matrix anchor and scenario-count staleness
# ---------------------------------------------------------------------

nAnchors <- 0L
nSymbols <- 0L
nCounts <- 0L

matrixRel <- "docs/design/feature-matrix.md"
matrixFile <- p(matrixRel)

if (!file.exists(matrixFile)) {
  report(sprintf("MATRIX: %s is missing", matrixRel))
} else {
  # Path aliases, mirroring the "Path aliases used in anchors" block at the
  # head of feature-matrix.md; keep in sync if that block changes.
  ALIAS_PATH <- c(
    RIB = "src/R_interface_bartcore.cpp",
    CAPI = "inst/include/dbarts/dbarts.h",
    MOD = "src/bartcore/model.hpp",
    CH = "src/bartcore/chain.hpp",
    FAC = "src/bartcore/facade.hpp",
    COM = "src/bartcore/combiner.hpp",
    MOV = "src/bartcore/moves.hpp",
    SAM = "src/bartcore/sampler.hpp"
  )
  R_ALIAS_FILES <- c(
    "bart.R",
    "dbarts.R",
    "spec.R",
    "rbart.R",
    "xbart.R",
    "data.R",
    "generics.R",
    "A_class.R",
    "bartcore.R"
  )
  RD_ALIAS <- c(
    sampler.Rd = "man/dbartsSampler-class.Rd",
    bart.Rd = "man/bart.Rd"
  )

  # Repo-tracked file inventory, for resolving anchors cited by bare
  # basename (e.g. "chain.hpp:1130", "bcf.md:370") rather than by alias.
  trackedFiles <- tryCatch(
    suppressWarnings(system2(
      "git",
      c("-C", root, "ls-files"),
      stdout = TRUE,
      stderr = FALSE
    )),
    error = function(e) character(0)
  )
  if (length(trackedFiles) == 0L) {
    trackedFiles <- list.files(root, recursive = TRUE, full.names = FALSE)
    trackedFiles <- trackedFiles[!grepl("(^|/)\\.git/", trackedFiles)]
  }
  basenameIndex <- split(trackedFiles, basename(trackedFiles))

  resolveToken <- function(tok) {
    if (tok %in% names(ALIAS_PATH)) {
      return(unname(ALIAS_PATH[tok]))
    }
    if (tok %in% names(RD_ALIAS)) {
      return(unname(RD_ALIAS[tok]))
    }
    if (tok %in% R_ALIAS_FILES) {
      return(file.path("R", tok))
    }
    if (tok == "TODO") {
      return("TODO")
    }
    if (grepl("/", tok, fixed = TRUE)) {
      return(tok)
    }
    hit <- basenameIndex[[tok]]
    if (is.null(hit)) {
      return(NA_character_)
    }
    hit[1L]
  }

  fileCache <- new.env(parent = emptyenv())
  getContent <- function(relPath) {
    key <- relPath
    if (!exists(key, envir = fileCache, inherits = FALSE)) {
      full <- p(relPath)
      val <- if (file.exists(full)) readLines(full, warn = FALSE) else NA
      assign(key, val, envir = fileCache)
    }
    get(key, envir = fileCache, inherits = FALSE)
  }
  fileExistsCached <- function(relPath) !identical(getContent(relPath), NA)

  # A cited symbol counts as present if the literal string occurs verbatim,
  # or - for a "Class::member" citation, the common case for an in-class C++
  # definition, which never repeats the class qualifier at its own
  # declaration - if both the class name and the member name occur
  # somewhere in the file (not necessarily adjacent). This trades a little
  # precision (it cannot tell "GaussianResponse::setActiveRows" from another
  # class's setActiveRows sharing the file) for not flagging every ordinary
  # in-class method as a dead symbol.
  symbolPresent <- function(sym, content) {
    if (any(grepl(sym, content, fixed = TRUE))) {
      return(TRUE)
    }
    parts <- strsplit(sym, "::", fixed = TRUE)[[1L]]
    if (length(parts) <= 1L) {
      return(FALSE)
    }
    all(vapply(
      parts,
      function(part) any(grepl(part, content, fixed = TRUE)),
      logical(1L)
    ))
  }

  matrixLines <- readLines(matrixFile, warn = FALSE)

  # --- 2a: ALIAS:line[-line] and path.R:line anchors, each optionally
  # preceded on the same line by a `symbol` naming what the line contains ---
  anchorRe <- paste0(
    "\\b((?:RIB|CAPI|MOD|CH|FAC|COM|MOV|SAM|TODO)|",
    "[A-Za-z_][A-Za-z0-9_./-]*\\.(?:R|Rd|md|hpp|cpp|c|h)):",
    "([0-9]+)(?:-([0-9]+))?"
  )

  for (i in seq_along(matrixLines)) {
    line <- matrixLines[i]
    m <- gregexpr(anchorRe, line, perl = TRUE)[[1]]
    if (m[1L] == -1L) {
      next
    }
    starts <- as.integer(m)
    lens <- attr(m, "match.length")
    capStart <- attr(m, "capture.start")
    capLen <- attr(m, "capture.length")
    for (k in seq_along(starts)) {
      tok <- substr(line, capStart[k, 1L], capStart[k, 1L] + capLen[k, 1L] - 1L)
      lo <- as.integer(substr(
        line,
        capStart[k, 2L],
        capStart[k, 2L] + capLen[k, 2L] - 1L
      ))
      hi <- if (capLen[k, 3L] > 0L) {
        as.integer(substr(
          line,
          capStart[k, 3L],
          capStart[k, 3L] + capLen[k, 3L] - 1L
        ))
      } else {
        lo
      }

      relPath <- resolveToken(tok)
      if (is.na(relPath)) {
        next
      } # outside this guard's alias/basename coverage
      nAnchors <- nAnchors + 1L
      content <- getContent(relPath)
      if (identical(content, NA)) {
        report(sprintf(
          "%s:%d: anchor '%s:%d' -> %s does not exist",
          matrixRel,
          i,
          tok,
          lo,
          relPath
        ))
        next
      }
      if (hi > length(content)) {
        report(sprintf(
          "%s:%d: anchor '%s:%s' -> %s has only %d lines",
          matrixRel,
          i,
          tok,
          if (hi != lo) paste0(lo, "-", hi) else as.character(lo),
          relPath,
          length(content)
        ))
      }

      # Skip the lookback when the anchor is itself backtick-wrapped
      # ("`spec.R:440`"): the tick right before it is that wrap's own
      # opener, not a preceding symbol span's closer, and pairing it with
      # an earlier tick would grab a run of prose instead of a symbol.
      matchEnd <- starts[k] + lens[k]
      selfWrapped <- substr(line, matchEnd, matchEnd) == "`" &&
        substr(line, starts[k] - 1L, starts[k] - 1L) == "`"
      if (!selfWrapped) {
        prefix <- sub("\\s+$", "", substr(line, 1L, starts[k] - 1L))
        symMatch <- regmatches(prefix, regexec("`([^`]+)`$", prefix))[[1L]]
        if (length(symMatch) == 2L) {
          sym <- symMatch[2L]
          nSymbols <- nSymbols + 1L
          if (!symbolPresent(sym, content)) {
            report(sprintf(
              "%s:%d: symbol `%s` (cited beside %s:%d) not found in %s",
              matrixRel,
              i,
              sym,
              tok,
              lo,
              relPath
            ))
          }
        }
      }
    }
  }

  # --- 2b: the rarer "path.R (symbol)" anchor with no line number ---
  symOnlyRe <- "\\b([A-Za-z0-9_./-]+\\.(?:R|Rd))\\s*\\(([A-Za-z_][A-Za-z0-9_:.]*)\\)"
  for (i in seq_along(matrixLines)) {
    line <- matrixLines[i]
    m <- gregexpr(symOnlyRe, line, perl = TRUE)[[1L]]
    if (m[1L] == -1L) {
      next
    }
    starts <- as.integer(m)
    capStart <- attr(m, "capture.start")
    capLen <- attr(m, "capture.length")
    for (k in seq_along(starts)) {
      tok <- substr(line, capStart[k, 1L], capStart[k, 1L] + capLen[k, 1L] - 1L)
      sym <- substr(line, capStart[k, 2L], capStart[k, 2L] + capLen[k, 2L] - 1L)
      relPath <- resolveToken(tok)
      if (is.na(relPath)) {
        next
      }
      content <- getContent(relPath)
      if (identical(content, NA)) {
        report(sprintf(
          "%s:%d: reference '%s' does not exist",
          matrixRel,
          i,
          relPath
        ))
        next
      }
      nSymbols <- nSymbols + 1L
      if (!symbolPresent(sym, content)) {
        report(sprintf(
          "%s:%d: symbol `%s` (cited beside %s) not found in %s",
          matrixRel,
          i,
          sym,
          tok,
          relPath
        ))
      }
    }
  }

  # --- 2c: recomputed scenario-count claims (e.g. footnote [f39]) ---
  # Each baseline .rds is recorded by a benchmarks/R/*.R harness that builds
  # its scenario list as a sequence of `result$<name> <- ...` assignments
  # inside one named function; the matrix cites how many scenarios that
  # baseline carries, and this recounts the names the harness itself defines.
  SCENARIO_SOURCES <- list(
    list(
      prefix = "multinomial-equivalence-",
      file = "benchmarks/R/multinomial-equivalence.R",
      fn = "runScenarios"
    ),
    list(
      prefix = "bcf-equivalence-",
      file = "benchmarks/R/bcf-equivalence.R",
      fn = "runScenarios"
    ),
    list(
      prefix = "equivalence-",
      file = "benchmarks/R/equivalence.R",
      fn = "makeScenarios"
    )
  )

  scenarioNamesIn <- function(relFile, fnName) {
    lines <- getContent(relFile)
    if (identical(lines, NA)) {
      return(NA_character_)
    }
    startLine <- grep(paste0("^", fnName, " *<- *function"), lines)
    if (length(startLine) == 0L) {
      return(NA_character_)
    }
    startLine <- startLine[1L]
    closes <- grep("^\\}", lines)
    closes <- closes[closes > startLine]
    endLine <- if (length(closes) > 0L) closes[1L] else length(lines)
    body <- lines[startLine:endLine]
    nm <- gregexpr("(?<=result\\$)[A-Za-z0-9_.]+(?=\\s*<-)", body, perl = TRUE)
    unique(unlist(regmatches(body, nm)))
  }

  blob <- paste(matrixLines, collapse = "\n")
  baselinePat <- "`([a-zA-Z0-9_.-]+\\.rds)`\\s*\\((\\d+)(?:\\s+scenarios)?\\)"
  matchedStrs <- regmatches(blob, gregexpr(baselinePat, blob, perl = TRUE))[[
    1L
  ]]
  for (s in matchedStrs) {
    g <- regmatches(s, regexec(baselinePat, s, perl = TRUE))[[1L]]
    fname <- g[2L]
    claimed <- as.integer(g[3L])
    src <- Find(function(e) startsWith(fname, e$prefix), SCENARIO_SOURCES)
    if (is.null(src)) {
      next
    } # not one of the recomputable baselines this guard knows
    nCounts <- nCounts + 1L

    baselineRel <- file.path("benchmarks/baselines", fname)
    if (!fileExistsCached(baselineRel)) {
      report(sprintf(
        "MATRIX: baseline '%s' cited (%d scenarios) does not exist at %s",
        fname,
        claimed,
        baselineRel
      ))
      next
    }
    names <- scenarioNamesIn(src$file, src$fn)
    if (identical(names, NA_character_)) {
      report(sprintf(
        "MATRIX: cannot locate %s() in %s to recount scenarios for '%s'",
        src$fn,
        src$file,
        fname
      ))
      next
    }
    actual <- length(names)
    if (actual != claimed) {
      report(sprintf(
        "MATRIX: '%s' claims %d scenarios but %s()'s %s has %d",
        fname,
        claimed,
        src$fn,
        src$file,
        actual
      ))
    }
  }
}

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

if (length(find) > 0L) {
  cat(paste0("FAIL: ", find), sep = "\n")
  cat(sprintf("check-doc-freshness: %d failure(s)\n", length(find)))
  quit(status = 1L, save = "no")
}
cat(sprintf(
  "check-doc-freshness: OK (%d design docs, %d plan docs indexed; %d anchors, %d symbols, %d scenario-count claims verified)\n",
  nDesign,
  nPlans,
  nAnchors,
  nSymbols,
  nCounts
))
quit(status = 0L, save = "no")
