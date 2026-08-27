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
# 3. AT-LINE ANCHOR CHECK: every file:line / file:line-line / file:line,
#    line anchor in EVERY docs/design/*.md (not just feature-matrix.md)
#    resolves against the live tree, and an identifier the doc cites
#    beside it must occur near the CITED line, not merely somewhere in
#    the file. EVERY backticked token in the anchor's own sentence (or
#    table cell) is a candidate, since prose pairs lists of names with
#    lists of anchors; notation (whitespace, brackets, parentheses,
#    template args, operators) is not a candidate at all; a candidate
#    also lands when it names the definition ENCLOSING the cited line.
#    STRICT (the nearest candidate immediately precedes the anchor, at
#    most three words between, and the target file defines that name
#    somewhere): a miss is a FAILURE - the name is in the file, just not
#    where the doc says. ADVISORY (the nearest candidate sits farther
#    back, is possessive, or is a name the target file never defines -
#    an R5 method named beside the bridge entry that refuses it, say):
#    a miss is only a WARNING. Bare filenames resolve under R/, src/,
#    inst/, tests/, benchmarks/, tools/, man/ (or docs/design,
#    docs/plans for a *.md target); ambiguous ones are skipped, counted.
#    A bare ":N" continuation inherits its own table cell's file, then
#    the file the table assigns that column, then the row's, then the
#    column's, then the prose run's. Doc-to-doc anchors are checked for
#    existence only. Skipped, not failed: fenced code blocks; a file or
#    section range a doc itself marks not-live/frozen/historical;
#    anchors annotated `retired:`/`unresolved` nearby; peer-package
#    citations (paths that do not resolve in this tree).
#
# 4. HASH RESOLVABILITY: every 7-12 character hex token that looks like a
#    commit reference in docs/design/*.md, TODO, and
#    benchmarks/baselines/MANIFEST must resolve via
#    `git cat-file -e <hash>^{commit}`. Excludes the hex part of a
#    "0x..." API hash literal, a hash tagged nearby as pre-rebase, a
#    hash whose own paragraph or table row names another repository
#    (stan4bart, bartCause, pymc-bart, the dbarts-1.0 compat branch and
#    the rest), and (locally only - CI's checkout is shallow, so this
#    whole check is skipped there) anything unresolvable simply because
#    history was not fetched.
#
# 5. INDEX STATUS: for every docs/design/*.md carrying a literal
#    "Status:" line in its first 12 lines, the leading status phrase
#    INDEX.md gives that doc must be confirmable, token by token,
#    case-insensitively, somewhere in the doc's own first-12-lines
#    block - catching either side going stale after a verdict changes.
#    Docs without a Status: line are counted and skipped.
#
# This is a DEAD-REFERENCE and DRIFTED-COUNT detector, not a semantic
# reviewer: it never re-adjudicates a SHIPPED/REFUSED/etc. cell VALUE, and
# Part 2's symbol check only asks whether the string still occurs in the
# file ANYWHERE, not at the cited line - a line number that has merely
# moved within a file that still exists and still contains the symbol is
# not flagged there (Part 3 covers that). An anchor this guard's
# alias/basename resolution cannot place is skipped, not failed.
#
# Usage: Rscript tools/check-doc-freshness.R [pkg-root]
# Exit status: 0 if every check passes, 1 if any fails (failures print
# with a leading "FAIL:" so they are easy to grep; ADVISORY misses print
# with a leading "WARN:" and do not affect the exit status).

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[1L] else "."
p <- function(...) file.path(root, ...)

find <- character(0L)
report <- function(msg) find <<- c(find, msg)

warnFind <- character(0L)
reportWarn <- function(msg) warnFind <<- c(warnFind, msg)

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
# Shared: path aliases, repo file inventory, and a per-file content
# cache. Used by Part 2 (feature-matrix only) and by Part 3 (every
# docs/design/*.md).
# ---------------------------------------------------------------------

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
  if (!is.null(hit) && file.exists(p(hit[1L]))) {
    return(hit[1L])
  }
  # A plan doc's tracked-index location can be stale relative to the
  # working tree (e.g. archived out from docs/plans/ into
  # docs/plans/archive/ after the index was snapshotted); try both
  # literal locations before giving up on the index's own answer.
  if (grepl("\\.md$", tok)) {
    planPath <- file.path("docs/plans", tok)
    if (file.exists(p(planPath))) {
      return(planPath)
    }
    archivePath <- file.path("docs/plans/archive", tok)
    if (file.exists(p(archivePath))) {
      return(archivePath)
    }
  }
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
# Part 3: at-line anchor check, across every docs/design/*.md
# ---------------------------------------------------------------------
#
# Unlike Part 2 (feature-matrix.md only, symbol-present-somewhere-in-file),
# this walks every design doc and, for each file:line / file:line-line /
# file:line,line,... anchor, requires an identifier the doc cites beside
# the anchor to occur within +/-2 lines of the CITED range, not merely
# somewhere in the file.
#
# CANDIDATES. Every backticked token in the anchor's own sentence (or
# table cell) is a candidate, not merely the nearest one: prose routinely
# pairs a slash- or comma-joined LIST of names with one anchor
# ("`refuseBartRedirectedFamily`/`bartRedirectedFamilies` (bart.R:2632)"),
# and only one member of the list can be the cited line's own name. The
# anchor passes if ANY candidate lands. Tokens that are notation rather
# than identifiers - anything carrying whitespace, brackets, parentheses
# (other than a trailing "()"), template arguments or operators, e.g.
# "[node.begin, node.end)", "muByTree[t][leafOf[t * n + i]]", "w_i * a^2",
# "family = \"probit\"" - are not candidates at all: they are prose, and no
# file contains them verbatim. A leading "Class::" / "obj$" / "obj@"
# qualifier and a trailing "()" are stripped before the search, since a
# definition never repeats them.
#
# ENCLOSING SCOPE. A candidate also lands when it names the definition
# ENCLOSING the cited line - "`Sampler` SAM:1583" cites a member of class
# Sampler, "`dbarts()` dbarts.R:449" a line inside that function body.
# The search walks outward from the cited line by indentation, so a
# sibling definition (whose block already closed) is never mistaken for an
# enclosing one.
#
# LIST TO LIST. When one sentence carries several anchors AND several
# candidates ("`bartcore_setData`, `bartcore_setModel`, RIB:4639/:4932"),
# the pairing is positional and a name can follow the anchor it belongs
# to; each anchor then passes if any candidate anywhere in the sentence
# lands.
#
# TIERS. STRICT (the nearest candidate immediately precedes the anchor, at
# most three words between): a miss is a FAILURE. ADVISORY (the nearest
# candidate sits farther back in the sentence): a miss is only a WARNING.
# An anchor with no candidate at all is not checked either way (nothing to
# compare). Anchors into another docs/design or docs/plans file are
# checked for existence/line-count only - identifier presence is
# meaningless across a doc-to-doc citation.
#
# A bare ":N"/"N-M"/"N,M" continuation (no filename) inherits, in order:
# the last file named in its own table CELL; the file the table assigns to
# that COLUMN (either a header-row cell that is itself a filename, or a
# prose line before the table that says which file each column's bare
# anchors belong to - see multinomial-mutation-arc.md's capability table,
# "the bridge column is `src/R_interface_bartcore.cpp`, the R-wrapper
# column `R/bartcore.R`"); the last file named earlier in the same ROW;
# the last file named in that column in an earlier row; and, failing all
# of those, the last file named in the prose run before the table. Outside
# a table it inherits the last file named earlier in the same run of prose
# (reset at each heading). A prose line that assigns files to columns
# names them ON BEHALF OF the table, so it does not itself become the
# prose run's "last named file".
#
# Skipped, not failed: anchors inside fenced code blocks; anchors in a
# file whose own head marks its citations as not-live (e.g. data-store.md,
# model-space-survey.md: "Code citations are at <hash>; they are not
# live."); anchors inside a doc-marked frozen/period-correct section
# range (e.g. multinomial-mutation-arc.md's "Sections 5-8 ... are not
# re-resynced here (they are not live pointers into the current tree)");
# anchors annotated `retired:` or `unresolved` nearby; bare filenames that
# do not resolve under R/, src/, inst/, tests/, benchmarks/, tools/, man/
# (or under docs/design, docs/plans for a *.md target) - this is how a
# peer-package citation (a path that simply is not in this tree) is told
# apart from a real dead anchor; and bare filenames that resolve to more
# than one tracked file (ambiguous - counted, not guessed at).

nAtLineAnchors <- 0L
nAtLineStrict <- 0L
nAtLineAdvisory <- 0L
nAmbiguousAnchor <- 0L

EXT_RE <- "(?:R|Rd|md|hpp|cpp|c|h)"
NAMED_ANCHOR_RE <- paste0(
  "\\b(RIB|CAPI|MOD|CH|FAC|COM|MOV|SAM|TODO|",
  "[A-Za-z_][A-Za-z0-9_./-]*\\.",
  EXT_RE,
  "):",
  "([0-9]+(?:[-,][0-9]+)*)"
)
# A bare ":N" continuation: a colon not immediately preceded by a
# filename-ish character (so it does not re-match the colon a named
# anchor already consumed) or by ")" (so a "29(2):405-417" journal
# volume(issue):pages citation is not mistaken for one).
BARE_ANCHOR_RE <- "(?<![A-Za-z0-9_./)]):([0-9]+(?:[-,][0-9]+)*)"
# A bare filename mention with no line number at all - just names the
# file the next bare ":N" continuations inherit.
NAMING_MENTION_RE <- paste0(
  "`(RIB|CAPI|MOD|CH|FAC|COM|MOV|SAM|TODO|",
  "[A-Za-z_][A-Za-z0-9_./-]*\\.",
  EXT_RE,
  ")`(?!:[0-9])"
)
# "the bridge column is `src/R_interface_bartcore.cpp`, the R-wrapper
# column `R/bartcore.R`" - a prose line assigning a file to a table
# column by that column's header label.
COL_ASSIGN_RE <- paste0(
  "(?:the\\s+)?([A-Za-z][A-Za-z0-9_.-]*(?:[ -][A-Za-z][A-Za-z0-9_.-]*)?)",
  "\\s+columns?\\s+(?:is\\s+|are\\s+)?`([^`]+)`"
)

BARE_DIRS_RE <- "^(R|src|inst|tests|benchmarks|tools|man)/"
DOC_DIRS_RE <- "^(docs/design|docs/plans)/"
bareBasenameIndex <- split(
  trackedFiles[grepl(BARE_DIRS_RE, trackedFiles)],
  basename(trackedFiles[grepl(BARE_DIRS_RE, trackedFiles)])
)
docBasenameIndex <- split(
  trackedFiles[grepl(DOC_DIRS_RE, trackedFiles)],
  basename(trackedFiles[grepl(DOC_DIRS_RE, trackedFiles)])
)

resolveAnchorTokenGeneral <- function(tok) {
  if (tok %in% names(ALIAS_PATH)) {
    return(unname(ALIAS_PATH[tok]))
  }
  if (tok %in% names(RD_ALIAS)) {
    return(unname(RD_ALIAS[tok]))
  }
  if (identical(tok, "TODO")) {
    return("TODO")
  }
  if (grepl("/", tok, fixed = TRUE)) {
    # A path with a directory component is only "resolved" if it is
    # actually in this tree; otherwise treat it as a peer-package
    # citation (e.g. a survey doc quoting another package's own
    # R/model.R) rather than reporting a false dead anchor. The
    # tracked-file index can itself be stale relative to the working
    # tree (e.g. an archive move the index was snapshotted before), so
    # a literal on-disk hit still counts as resolved.
    if (tok %in% trackedFiles || file.exists(p(tok))) {
      return(tok)
    }
    return(NA_character_)
  }
  ext <- tolower(sub(".*\\.", "", tok))
  hits <- if (identical(ext, "md")) {
    docBasenameIndex[[tok]]
  } else {
    bareBasenameIndex[[tok]]
  }
  if (!is.null(hits)) {
    uniqueHits <- unique(hits)
    if (length(uniqueHits) > 1L) {
      nAmbiguousAnchor <<- nAmbiguousAnchor + 1L
      return(NA_character_)
    }
    if (file.exists(p(uniqueHits[1L]))) {
      return(uniqueHits[1L])
    }
  }
  # As in resolveToken(): a plan doc's tracked-index location can be
  # stale relative to the working tree after an archive move.
  if (identical(ext, "md")) {
    planPath <- file.path("docs/plans", tok)
    if (file.exists(p(planPath))) {
      return(planPath)
    }
    archivePath <- file.path("docs/plans/archive", tok)
    if (file.exists(p(archivePath))) {
      return(archivePath)
    }
  }
  if (is.null(hits)) {
    return(NA_character_)
  }
  unique(hits)[1L]
}

countWords <- function(s) {
  mm <- gregexpr("[A-Za-z0-9_]+", s, perl = TRUE)[[1L]]
  if (mm[1L] == -1L) 0L else length(mm)
}

parseLineSpec <- function(spec) {
  parts <- strsplit(spec, ",", fixed = TRUE)[[1L]]
  los <- integer(length(parts))
  his <- integer(length(parts))
  for (j in seq_along(parts)) {
    if (grepl("-", parts[j], fixed = TRUE)) {
      rp <- as.integer(strsplit(parts[j], "-", fixed = TRUE)[[1L]])
      los[j] <- min(rp)
      his[j] <- max(rp)
    } else {
      v <- as.integer(parts[j])
      los[j] <- v
      his[j] <- v
    }
  }
  list(lo = min(los), hi = max(his))
}

# --- candidate identifiers -------------------------------------------

# An anchor is not an identifier to search for: a bare list of citations
# ("`R/data.R:229-318`, `R/model.R:93-133`") carries no symbol at all.
ANCHOR_TOKEN_RE <- "^\\(?[A-Za-z0-9_./-]*:[0-9]+(?:[-,][0-9]+)*\\)?$"
# Nor is a bare filename ("applied across `R/data.R`, `R/A_class.R`").
FILENAME_TOKEN_RE <- paste0("^[A-Za-z0-9_./-]+\\.", EXT_RE, "$")
# Notation, not an identifier: whitespace, brackets, parentheses,
# template args, operators, quotes, path separators.
NOTATION_RE <- "[\\s<>\\[\\](){}+*,=|\"'^!&%;?~/\\\\]"
IDENT_SHAPE_RE <- "^[A-Za-z_][A-Za-z0-9_.-]*$"

normalizeIdent <- function(tok) {
  t <- trimws(tok)
  t <- sub("\\(\\)$", "", t)
  repeat {
    t2 <- sub("^[A-Za-z_.][A-Za-z0-9_.]*(?:::|\\$|@)", "", t, perl = TRUE)
    t2 <- sub("^[$@]", "", t2, perl = TRUE)
    if (identical(t2, t)) {
      break
    }
    t <- t2
  }
  t
}

# The identifier a backticked token contributes, or NA when the token is
# an anchor, a filename, or notation rather than a name.
identOf <- function(tok) {
  if (grepl(ANCHOR_TOKEN_RE, trimws(tok), perl = TRUE)) {
    return(NA_character_)
  }
  if (grepl(FILENAME_TOKEN_RE, trimws(tok), perl = TRUE)) {
    return(NA_character_)
  }
  ident <- normalizeIdent(tok)
  if (grepl(NOTATION_RE, ident, perl = TRUE)) {
    return(NA_character_)
  }
  if (!grepl(IDENT_SHAPE_RE, ident, perl = TRUE) || nchar(ident) < 2L) {
    return(NA_character_)
  }
  ident
}

# Every backticked token in `text`, as (ident, gap-to-end-of-text) pairs.
# When the anchor is itself the thing inside the nearest backtick pair
# ("`spec.R:440`"), that pair's own opener is dropped from the search so
# it is not paired with an earlier, unrelated closer.
collectCandidates <- function(text, selfWrapped = FALSE) {
  search <- text
  if (
    selfWrapped &&
      nchar(search) > 0L &&
      substr(search, nchar(search), nchar(search)) == "`"
  ) {
    search <- substr(search, 1L, nchar(search) - 1L)
  }
  m <- gregexpr("`([^`]+)`", search, perl = TRUE)[[1L]]
  if (m[1L] == -1L) {
    return(list())
  }
  lens <- attr(m, "match.length")
  out <- list()
  for (k in seq_along(m)) {
    full <- substr(search, m[k], m[k] + lens[k] - 1L)
    tok <- substr(full, 2L, nchar(full) - 1L)
    ident <- identOf(tok)
    if (is.na(ident)) {
      next
    }
    gap <- substr(text, m[k] + lens[k], nchar(text))
    out[[length(out) + 1L]] <- list(
      tok = tok,
      ident = ident,
      gapWords = countWords(gap),
      # "`dbartsSpec()`'s BCF composition (spec.R:658)", "`run`'s own
      # `:377` arms": a possessive names the OWNER of what the anchor
      # points at, not the name the cited line itself carries.
      possessive = grepl("^\\s*['\u2019]s?\\b", gap, perl = TRUE)
    )
  }
  out
}

# --- enclosing-scope acceptance --------------------------------------
#
# Per line: the name that line DEFINES (NA if none), and an effective
# indentation. A line that only closes a block counts as half a level
# further out than the block it closes, so scanning outward from the
# cited line never mistakes an already-closed sibling for an enclosing
# scope.

R_DEF_RE <- "^\\s*([A-Za-z_.][A-Za-z0-9_.]*)\\s*(?:<-|=)\\s*function\\b"
# An R assignment whose right-hand side opens and does not close on the
# same line ("anchorSamplers <- list(") encloses the lines that follow it
# exactly as a function definition does.
R_BLOCK_RE <- "^\\s*([A-Za-z_.][A-Za-z0-9_.]*)\\s*(?:<-|=)\\s*(?:[A-Za-z_.][A-Za-z0-9_.:]*)?\\(\\s*$"
CXX_TYPE_RE <- "^\\s*(?:template\\s*<[^>]*>\\s*)?(?:class|struct|union|enum(?:\\s+class)?)\\s+([A-Za-z_][A-Za-z0-9_]*)\\b"
CXX_FN_RE <- "^\\s*(?:[A-Za-z_][A-Za-z0-9_:<>,*&\\s]*[\\s*&])([A-Za-z_][A-Za-z0-9_]*)\\s*\\("
CXX_FN_KEYWORD_RE <- "\\b(?:return|if|for|while|switch|else|do|new|delete|throw|case|sizeof|catch)\\s*$"
DEFINE_RE <- "^\\s*#\\s*define\\s+([A-Za-z_][A-Za-z0-9_]*)"
CLOSER_RE <- "^\\s*[)}\\]]+[;,]?\\s*$"

scopeCache <- new.env(parent = emptyenv())
scopeInfo <- function(relPath) {
  if (!exists(relPath, envir = scopeCache, inherits = FALSE)) {
    content <- getContent(relPath)
    if (identical(content, NA)) {
      val <- NULL
    } else {
      def <- rep(NA_character_, length(content))
      hit <- grepl(R_DEF_RE, content, perl = TRUE)
      def[hit] <- sub(paste0(R_DEF_RE, ".*$"), "\\1", content[hit], perl = TRUE)
      hit <- is.na(def) & grepl(R_BLOCK_RE, content, perl = TRUE)
      def[hit] <- sub(
        paste0(R_BLOCK_RE, ".*$"),
        "\\1",
        content[hit],
        perl = TRUE
      )
      hit <- is.na(def) & grepl(CXX_TYPE_RE, content, perl = TRUE)
      def[hit] <- sub(
        paste0(CXX_TYPE_RE, ".*$"),
        "\\1",
        content[hit],
        perl = TRUE
      )
      hit <- is.na(def) & grepl(DEFINE_RE, content, perl = TRUE)
      def[hit] <- sub(
        paste0(DEFINE_RE, ".*$"),
        "\\1",
        content[hit],
        perl = TRUE
      )
      hit <- is.na(def) & grepl(CXX_FN_RE, content, perl = TRUE)
      if (any(hit)) {
        nm <- sub(paste0(CXX_FN_RE, ".*$"), "\\1", content[hit], perl = TRUE)
        # "return foo(x);" defines nothing.
        pre <- sub(
          "^(.*?)\\b[A-Za-z_][A-Za-z0-9_]*\\s*\\(.*$",
          "\\1",
          content[hit],
          perl = TRUE
        )
        nm[grepl(CXX_FN_KEYWORD_RE, pre, perl = TRUE)] <- NA_character_
        def[hit] <- nm
      }
      # Only a definition head or a block-closing line delimits a scope;
      # every other line (statement, comment, continuation, access
      # specifier, macro body) is transparent to the outward walk, so a
      # comment or a "public:" at column 0 cannot hide the class that
      # encloses the cited line.
      closer <- grepl(CLOSER_RE, content, perl = TRUE)
      ind <- rep(Inf, length(content))
      structural <- !is.na(def) | closer
      ind[structural] <- nchar(sub("^([ \t]*).*$", "\\1", content[structural]))
      ind[closer] <- ind[closer] - 0.5
      # The name a line LEADS with, past any list/comment marker: how a
      # plain-text file (TODO's "multinomial-counts-mutation:" entries)
      # introduces a name it has no syntax to define.
      head <- rep(NA_character_, length(content))
      leads <- grepl("^[\\s>*#/|+-]*[A-Za-z_.]", content, perl = TRUE)
      head[leads] <- sub(
        "^[\\s>*#/|+-]*([A-Za-z_.][A-Za-z0-9_.-]*).*$",
        "\\1",
        content[leads],
        perl = TRUE
      )
      val <- list(
        indent = ind,
        def = def,
        defined = unique(c(def[!is.na(def)], head[!is.na(head)]))
      )
    }
    assign(relPath, val, envir = scopeCache)
  }
  get(relPath, envir = scopeCache, inherits = FALSE)
}

# Does the file define this name anywhere at all? A name the target file
# never defines (a capability named on another surface - "`setWeights`
# (RIB:2661-2665)" names the R5 method the bridge REFUSES there, not a C
# symbol; a C++ predicate named beside the tinytest range that pins it)
# cannot be located by line at all, so its absence from the cited window
# is not evidence a line number drifted. Those miss as ADVISORY; a
# STRICT failure means the name IS defined in the file, just not where
# the doc says.
definedInFile <- function(relPath, ident) {
  si <- scopeInfo(relPath)
  !is.null(si) && ident %in% si$defined
}

# Does some definition ENCLOSING the cited line carry this name? Walks
# outward: a line qualifies as an enclosing construct only when its
# effective indentation is strictly less than every line between it and
# the cited one.
enclosingDefines <- function(relPath, lineNo, ident) {
  si <- scopeInfo(relPath)
  if (is.null(si) || lineNo < 1L || lineNo > length(si$indent)) {
    return(FALSE)
  }
  ind <- rev(si$indent[seq_len(lineNo)])
  if (length(ind) < 2L) {
    return(FALSE)
  }
  outer <- ind[-1L] < cummin(ind)[-length(ind)]
  if (!any(outer)) {
    return(FALSE)
  }
  at <- lineNo - which(outer)
  any(!is.na(si$def[at]) & si$def[at] == ident)
}

# --- sentence scoping -------------------------------------------------
#
# A sentence ends at ".", "!" or "?" followed by whitespace; the common
# abbreviations and every "file.R:12" style anchor are excluded (an
# extension is followed by a colon or a letter, never by a space).
SENT_END_RE <- "(?<!\\be\\.g)(?<!\\bi\\.e)(?<!\\bvs)(?<!\\bcf)(?<!\\betc)(?<!\\bal)(?<!\\bsec)(?<!\\beq)[.!?](?=\\s)"

# Every file:line / bare :line anchor inside a scoped span of text.
anchorsIn <- function(text) {
  hits <- integer(0)
  for (re in c(NAMED_ANCHOR_RE, BARE_ANCHOR_RE)) {
    m <- gregexpr(re, text, perl = TRUE)[[1L]]
    if (m[1L] != -1L) {
      hits <- c(hits, as.integer(m))
    }
  }
  hits
}

sentenceEnds <- function(text) {
  m <- gregexpr(SENT_END_RE, text, perl = TRUE)[[1L]]
  if (m[1L] == -1L) integer(0) else as.integer(m)
}

designFilesRel <- if (dir.exists(p("docs/design"))) {
  sort(list.files(p("docs/design"), pattern = "\\.md$"))
} else {
  character(0)
}

for (df in designFilesRel) {
  relDoc <- file.path("docs/design", df)
  lines <- getContent(relDoc)
  if (identical(lines, NA) || length(lines) == 0L) {
    next
  }
  n <- length(lines)

  headSpan <- lines[seq_len(min(10L, n))]
  if (any(grepl("not live", headSpan, ignore.case = TRUE))) {
    next
  }

  # Frozen (period-correct / not-live) section ranges, as the doc itself
  # marks them.
  blob <- paste(lines, collapse = " ")

  # The bare word "plan" before a bare ":N" ("plan :4852-4854", "landed
  # per :4923-4924" continuing the same thought) is shorthand for this
  # doc's own companion docs/plans/*.md file, almost always named
  # in-text near the top (e.g. "(docs/plans/multiforest-extension-
  # surface.md, M4)") without backticks around it, so the ordinary
  # naming-mention scan (which requires them) never sees it.
  planDocForDoc <- {
    pm <- regmatches(
      blob,
      regexpr("docs/plans/(?:archive/)?[A-Za-z0-9_.+-]+\\.md", blob, perl = TRUE)
    )
    if (length(pm) == 1L && nzchar(pm) && file.exists(p(pm))) {
      pm
    } else {
      NA_character_
    }
  }

  frozen <- list()
  fm <- gregexpr("[Ss]ections? +([0-9]+)-([0-9]+)", blob, perl = TRUE)[[1L]]
  if (fm[1L] != -1L) {
    caps <- attr(fm, "capture.start")
    capLens <- attr(fm, "capture.length")
    mlens <- attr(fm, "match.length")
    headingLineNo <- grep("^#+\\s", lines)
    headingNumLineNo <- grep("^#+\\s+[0-9]+\\.", lines)
    headingNumVal <- as.integer(sub(
      "^#+\\s+([0-9]+)\\..*$",
      "\\1",
      lines[headingNumLineNo]
    ))
    for (k in seq_along(fm)) {
      loNum <- as.integer(substr(
        blob,
        caps[k, 1L],
        caps[k, 1L] + capLens[k, 1L] - 1L
      ))
      hiNum <- as.integer(substr(
        blob,
        caps[k, 2L],
        caps[k, 2L] + capLens[k, 2L] - 1L
      ))
      after <- substr(
        blob,
        fm[k] + mlens[k],
        min(nchar(blob), fm[k] + mlens[k] + 300L)
      )
      if (
        !grepl(
          "not live|period-correct|not re-resynced|not resynced",
          after,
          ignore.case = TRUE
        )
      ) {
        next
      }
      startIdx <- headingNumLineNo[headingNumVal == loNum]
      if (length(startIdx) == 0L) {
        next
      }
      hiIdx <- headingNumLineNo[headingNumVal == hiNum]
      hiLine <- if (length(hiIdx) > 0L) hiIdx[1L] else startIdx[1L]
      afterHeadings <- headingLineNo[headingLineNo > hiLine]
      endLine <- if (length(afterHeadings) > 0L) afterHeadings[1L] - 1L else n
      frozen[[length(frozen) + 1L]] <- c(startIdx[1L], endLine)
    }
  }
  inFrozen <- function(ln) {
    if (length(frozen) == 0L) {
      return(FALSE)
    }
    any(vapply(frozen, function(r) ln >= r[1L] && ln <= r[2L], logical(1L)))
  }

  lastFile <- NA_character_
  # Tracks only the last-named CODE/source file (never a docs/design or
  # docs/plans target) - what a freshly-entered table's per-column
  # fallback seeds from. Flowing prose right before a table often names
  # a doc alongside the code file the table is actually about (e.g.
  # "against `src/foo.cpp` and `docs/design/bar.md:12-34`:"), and a
  # plain single "last named file" tracker would hand the table the doc
  # by mistake.
  lastCodeFile <- NA_character_
  inTable <- FALSE
  colExplicit <- list()
  colImplicit <- list()
  proseColAssign <- list()
  inFence <- FALSE

  for (i in seq_len(n)) {
    line <- lines[i]

    if (grepl("^#+\\s", line)) {
      lastFile <- NA_character_
      lastCodeFile <- NA_character_
      inTable <- FALSE
      colExplicit <- list()
      colImplicit <- list()
      proseColAssign <- list()
    }
    if (grepl("^\\s*```", line)) {
      inFence <- !inFence
      next
    }
    if (inFence || inFrozen(i)) {
      next
    }

    isTableLine <- grepl("^\\s*\\|.*\\|\\s*$", line)

    # A prose line that assigns files to a table's columns names those
    # files on the table's behalf: record the assignment, and keep the
    # files out of the plain-prose "last named file" tracking.
    assignedHere <- character(0)
    if (!isTableLine) {
      am <- gregexpr(COL_ASSIGN_RE, line, perl = TRUE)[[1L]]
      if (am[1L] != -1L) {
        caps <- attr(am, "capture.start")
        capLens <- attr(am, "capture.length")
        for (k in seq_along(am)) {
          label <- substr(line, caps[k, 1L], caps[k, 1L] + capLens[k, 1L] - 1L)
          ftok <- substr(line, caps[k, 2L], caps[k, 2L] + capLens[k, 2L] - 1L)
          ftok <- sub(":[0-9].*$", "", ftok)
          rpA <- resolveAnchorTokenGeneral(ftok)
          if (is.na(rpA) || grepl(DOC_DIRS_RE, rpA)) {
            next
          }
          key <- gsub("[^a-z0-9]", "", tolower(label))
          if (nzchar(key)) {
            proseColAssign[[key]] <- rpA
            assignedHere <- c(assignedHere, rpA)
          }
        }
      }
    }

    if (isTableLine && !inTable) {
      colExplicit <- list()
      colImplicit <- list()
      # A markdown table header can name a whole column by its bare
      # filename with no backticks and no line number at all (e.g.
      # "| generic | R/generics.R | default expression |"), or a prose
      # line just above can assign one to the column by its header
      # label ("the bridge column is `src/R_interface_bartcore.cpp`").
      cells <- strsplit(line, "|", fixed = TRUE)[[1L]]
      for (ci in seq_along(cells)) {
        if (ci == 1L) {
          next
        }
        cellTrim <- trimws(gsub("`", "", cells[ci], fixed = TRUE))
        rpHead <- resolveAnchorTokenGeneral(cellTrim)
        if (!is.na(rpHead) && !grepl(DOC_DIRS_RE, rpHead)) {
          colExplicit[[as.character(ci - 1L)]] <- rpHead
          # Also seeds the plain-prose fallback, so a paragraph right
          # after the table ("`predictBlend` (:727-738) joins the
          # forwarding sites") that keeps citing the same file the
          # header just named does not fall back to whatever code file
          # was last named BEFORE the table instead.
          lastFile <- rpHead
          lastCodeFile <- rpHead
          next
        }
        if (length(proseColAssign) > 0L) {
          key <- gsub("[^a-z0-9]", "", tolower(cellTrim))
          if (nzchar(key)) {
            for (lab in names(proseColAssign)) {
              if (
                identical(lab, key) ||
                  grepl(lab, key, fixed = TRUE) ||
                  grepl(key, lab, fixed = TRUE)
              ) {
                colExplicit[[as.character(ci - 1L)]] <- proseColAssign[[lab]]
                break
              }
            }
          }
        }
      }
    }
    inTable <- isTableLine

    events <- list()

    mm <- gregexpr(NAMING_MENTION_RE, line, perl = TRUE)[[1L]]
    if (mm[1L] != -1L) {
      caps <- attr(mm, "capture.start")
      capLens <- attr(mm, "capture.length")
      for (k in seq_along(mm)) {
        tok <- substr(line, caps[k, 1L], caps[k, 1L] + capLens[k, 1L] - 1L)
        events[[length(events) + 1L]] <- list(
          pos = mm[k],
          type = "name",
          tok = tok
        )
      }
    }

    nm <- gregexpr(NAMED_ANCHOR_RE, line, perl = TRUE)[[1L]]
    if (nm[1L] != -1L) {
      caps <- attr(nm, "capture.start")
      capLens <- attr(nm, "capture.length")
      mlens <- attr(nm, "match.length")
      for (k in seq_along(nm)) {
        tok <- substr(line, caps[k, 1L], caps[k, 1L] + capLens[k, 1L] - 1L)
        spec <- substr(line, caps[k, 2L], caps[k, 2L] + capLens[k, 2L] - 1L)
        events[[length(events) + 1L]] <- list(
          pos = nm[k],
          end = nm[k] + mlens[k] - 1L,
          type = "named",
          tok = tok,
          spec = spec
        )
      }
    }

    bm <- gregexpr(BARE_ANCHOR_RE, line, perl = TRUE)[[1L]]
    if (bm[1L] != -1L) {
      caps <- attr(bm, "capture.start")
      capLens <- attr(bm, "capture.length")
      mlens <- attr(bm, "match.length")
      for (k in seq_along(bm)) {
        spec <- substr(line, caps[k, 1L], caps[k, 1L] + capLens[k, 1L] - 1L)
        events[[length(events) + 1L]] <- list(
          pos = bm[k],
          end = bm[k] + mlens[k] - 1L,
          type = "bare",
          spec = spec
        )
      }
    }

    if (length(events) == 0L) {
      next
    }
    events <- events[order(vapply(events, function(e) e$pos, numeric(1L)))]

    # Cell boundaries (table) and sentence boundaries (prose), so an
    # anchor's candidate identifiers are scoped to the cell or sentence
    # it actually sits in.
    bars <- if (isTableLine) {
      as.integer(gregexpr("|", line, fixed = TRUE)[[1L]])
    } else {
      integer(0)
    }
    sents <- if (!isTableLine) sentenceEnds(line) else integer(0)

    # A line naming two files for two different purposes (e.g. "against
    # `src/foo.cpp` and `docs/design/bar.md:12-34`:", right before a
    # table about foo.cpp) must not let the doc mention win lastFile
    # just because it comes second - once this line has named a
    # code/source file, a later doc mention on the SAME line no longer
    # updates the plain-prose fallback (it is still fully checked on
    # its own terms).
    codeSetThisLine <- FALSE
    # In a table, a named file is scoped to its own CELL first, then
    # offered to the rest of the ROW; the column keeps its own memory.
    cellFile <- NA_character_
    cellKey <- ""
    rowFile <- NA_character_

    for (ev in events) {
      col <- if (isTableLine) {
        sum(bars < ev$pos)
      } else {
        NA_integer_
      }
      colKey <- as.character(col)
      if (isTableLine && !identical(colKey, cellKey)) {
        cellKey <- colKey
        cellFile <- NA_character_
      }

      if (ev$type %in% c("name", "named")) {
        rp <- resolveAnchorTokenGeneral(ev$tok)
        if (is.na(rp)) {
          # An explicit filename this pass cannot resolve in-tree (a
          # peer-package citation, e.g. "stan4bart's own loop
          # (`src/init.cpp:642`)") is strong evidence the surrounding
          # text has moved on to a file this repo does not have; a
          # later bare ":N" here must not silently fall back to some
          # unrelated file named much earlier instead.
          if (isTableLine) {
            cellFile <- NA_character_
            colImplicit[[colKey]] <- NA_character_
          } else {
            lastFile <- NA_character_
            lastCodeFile <- NA_character_
          }
          next
        }
        isDocTarget <- grepl(DOC_DIRS_RE, rp)
        if (isDocTarget) {
          # A doc-to-doc mention is scoped to its own cell in a table
          # (the capability table's "status" column carries one-off doc
          # citations that must not poison every later bare ":N" in
          # that column) and, in prose, seeds the fallback only when
          # this line has not already named a code file.
          if (isTableLine) {
            cellFile <- rp
          } else if (!codeSetThisLine) {
            lastFile <- rp
          }
        } else if (isTableLine) {
          cellFile <- rp
          rowFile <- rp
          colImplicit[[colKey]] <- rp
        } else if (!(rp %in% assignedHere)) {
          lastFile <- rp
          lastCodeFile <- rp
          codeSetThisLine <- TRUE
        }
        if (ev$type == "name") {
          next
        }
      } else {
        # A bare ":N"/"N-M" preceded by the bare word "plan"/"plans" is
        # a different, unmarked shorthand ("plan :4852-4854" = "see the
        # companion docs/plans/*.md file at that line") for this doc's
        # own companion plan doc, not a continuation of whatever file
        # was last named for an unrelated reason; resolve it there (and
        # let it seed lastFile like any other doc-target naming, so a
        # nearby continuation - "landed per :4923-4924" - inherits it
        # too), or skip if this doc names no docs/plans/*.md at all.
        precedingWord <- substr(line, max(1L, ev$pos - 8L), ev$pos - 1L)
        if (ev$pos <= 8L && i > 1L) {
          # Heavily-wrapped prose can break "plan" onto the previous
          # physical line from the ":N" it governs.
          precedingWord <- paste0(
            substr(
              lines[i - 1L],
              max(1L, nchar(lines[i - 1L]) - 8L),
              nchar(lines[i - 1L])
            ),
            precedingWord
          )
        }
        isPlanRef <- grepl("\\bplans?\\s*$", precedingWord, ignore.case = TRUE)
        if (isPlanRef) {
          if (is.na(planDocForDoc)) {
            next
          }
          # Unlike an incidental doc mention, "plan" is an explicit,
          # unambiguous marker the author chose deliberately - it wins
          # even over a code file named earlier on the same line (e.g.
          # "combiner.hpp:1096-1099; plan :1958-1978, landed :4928-4932"
          # names two DIFFERENT, independent things in sequence, not an
          # aside).
          if (!isTableLine) {
            lastFile <- planDocForDoc
          }
        }
        rp <- if (isPlanRef) {
          planDocForDoc
        } else if (!isTableLine) {
          lastFile
        } else if (!is.na(cellFile)) {
          cellFile
        } else if (!is.null(colExplicit[[colKey]])) {
          colExplicit[[colKey]]
        } else if (!is.na(rowFile)) {
          rowFile
        } else if (
          !is.null(colImplicit[[colKey]]) && !is.na(colImplicit[[colKey]])
        ) {
          colImplicit[[colKey]]
        } else {
          lastCodeFile
        }
        if (is.na(rp)) {
          next
        }
        if (isTableLine) {
          cellFile <- rp
        }
      }

      ctxLo <- max(1L, ev$pos - 300L)
      ctxHi <- min(nchar(line), ev$end + 300L)
      if (
        grepl(
          "\\bretired\\b|unresolved|historical|stale at",
          substr(line, ctxLo, ctxHi),
          ignore.case = TRUE
        )
      ) {
        next
      }

      nAtLineAnchors <- nAtLineAnchors + 1L
      content <- getContent(rp)
      if (identical(content, NA)) {
        report(sprintf(
          "%s:%d: anchor '%s' -> %s does not exist",
          relDoc,
          i,
          ev$spec,
          rp
        ))
        next
      }
      rng <- parseLineSpec(ev$spec)
      if (rng$hi > length(content)) {
        report(sprintf(
          "%s:%d: anchor '%s' -> %s has only %d lines",
          relDoc,
          i,
          ev$spec,
          rp,
          length(content)
        ))
      }

      if (grepl(DOC_DIRS_RE, rp)) {
        next # doc-to-doc: existence/line-count only, no identifier check
      }

      selfWrapped <- ev$pos > 1L &&
        substr(line, ev$pos - 1L, ev$pos - 1L) == "`" &&
        ev$end < nchar(line) &&
        substr(line, ev$end + 1L, ev$end + 1L) == "`"

      afterStart <- if (selfWrapped) ev$end + 2L else ev$end + 1L
      if (isTableLine) {
        cellStart <- if (any(bars < ev$pos)) {
          max(bars[bars < ev$pos]) + 1L
        } else {
          1L
        }
        cellEnd <- if (any(bars > ev$end)) {
          min(bars[bars > ev$end]) - 1L
        } else {
          nchar(line)
        }
        window <- substr(line, cellStart, ev$pos - 1L)
        after <- substr(line, afterStart, cellEnd)
      } else {
        before <- sents[sents < ev$pos]
        sentStart <- if (length(before) > 0L) max(before) + 1L else 1L
        afterEnds <- sents[sents > ev$end]
        sentEnd <- if (length(afterEnds) > 0L) min(afterEnds) else nchar(line)
        window <- substr(line, sentStart, ev$pos - 1L)
        after <- substr(line, afterStart, sentEnd)
        # A sentence that began on an earlier line: keep walking back
        # while this line carries no sentence break before the anchor.
        j <- i
        while (sentStart == 1L && j > 1L && (i - j) < 3L) {
          prevLine <- lines[j - 1L]
          if (
            grepl("^\\s*$", prevLine) ||
              grepl("^#+\\s", prevLine) ||
              grepl("^\\s*\\|", prevLine) ||
              grepl("^\\s*```", prevLine)
          ) {
            break
          }
          j <- j - 1L
          pe <- sentenceEnds(prevLine)
          if (length(pe) > 0L) {
            window <- paste0(
              substr(prevLine, max(pe) + 1L, nchar(prevLine)),
              " ",
              window
            )
            break
          }
          window <- paste0(prevLine, " ", window)
        }
      }

      if (rng$lo > length(content)) {
        next # already reported above as out-of-range; nothing sane to compare
      }

      # How many anchors this one sentence (or cell) carries, its own
      # included - the sentence can span lines, so this is counted over
      # the scoped text rather than over the physical line's events.
      nSentenceAnchors <- 1L +
        length(anchorsIn(window)) +
        length(anchorsIn(after))

      cands <- collectCandidates(window, selfWrapped)
      if (length(cands) == 0L && nSentenceAnchors < 2L) {
        next # no candidate identifier to compare against - not checked
      }

      lo2 <- max(1L, rng$lo - 2L)
      hi2 <- min(length(content), rng$hi + 2L)
      windowText <- gsub("\\s+", " ", paste(content[lo2:hi2], collapse = " "))
      lands <- function(cs) {
        for (cand in cs) {
          if (grepl(cand$ident, windowText, fixed = TRUE)) {
            return(TRUE)
          }
          if (enclosingDefines(rp, rng$lo, cand$ident)) {
            return(TRUE)
          }
        }
        FALSE
      }

      present <- lands(cands)
      # List to list: several anchors and several names in one sentence
      # pair off positionally, and a name can follow its own anchor.
      if (!present && nSentenceAnchors > 1L) {
        tailCands <- collectCandidates(after)
        if (length(cands) + length(tailCands) > 1L) {
          present <- lands(tailCands)
        }
      }
      if (length(cands) == 0L) {
        next
      }

      nearest <- cands[[length(cands)]]
      strict <- nearest$gapWords <= 3L &&
        !nearest$possessive &&
        (present || definedInFile(rp, nearest$ident))
      if (strict) {
        nAtLineStrict <- nAtLineStrict + 1L
        if (!present) {
          report(sprintf(
            "%s:%d: identifier `%s` (cited beside '%s') not found within %s:%d-%d",
            relDoc,
            i,
            nearest$tok,
            ev$spec,
            rp,
            lo2,
            hi2
          ))
        }
      } else {
        nAtLineAdvisory <- nAtLineAdvisory + 1L
        if (!present) {
          reportWarn(sprintf(
            "%s:%d: identifier `%s` (cited beside '%s') not found within %s:%d-%d",
            relDoc,
            i,
            nearest$tok,
            ev$spec,
            rp,
            lo2,
            hi2
          ))
        }
      }
    }
  }
}

# ---------------------------------------------------------------------
# Part 4: commit-hash resolvability
# ---------------------------------------------------------------------
#
# Every 7-12 character hex token that looks like a commit reference in
# docs/design/*.md, TODO, and benchmarks/baselines/MANIFEST must resolve
# with `git cat-file -e <hash>^{commit}`. A token must contain at least
# one a-f letter (otherwise it is far more likely a plain count/line
# number than a hash) and is excluded when it is really the hex part of
# a "0x..." API hash literal (an "x" is never a valid hex digit, so a
# 0x-prefixed run never has a \b boundary a 7-12 char submatch could
# start from) or is tagged nearby as a pre-rebase-candidate commit
# (the bartcore-pre-cran-rebase tag; such an object can legitimately be
# unreachable/pruned later). A hash whose own paragraph (or table row)
# names ANOTHER repository - stan4bart, bartCause, treatSens, bairrtt,
# pymc-bart/bartrs, bart-playground, bart-comp-efficiency, or the
# dbarts-1.0 compat branch - is a citation into that repository's
# history, which this clone cannot resolve and should not try to.
# Because CI checks out shallow, the whole check is skipped there - a
# shallow clone cannot resolve any hash.

nHashChecked <- 0L
isShallow <- tryCatch(
  identical(
    trimws(suppressWarnings(system2(
      "git",
      c("-C", root, "rev-parse", "--is-shallow-repository"),
      stdout = TRUE,
      stderr = FALSE
    ))),
    "true"
  ),
  error = function(e) FALSE
)

if (isShallow) {
  cat("hash check skipped (shallow clone)\n")
} else {
  hashRel <- c(
    file.path("docs/design", designFilesRel),
    "TODO",
    "benchmarks/baselines/MANIFEST"
  )
  hashRe <- "\\b[0-9a-fA-F]{7,12}\\b"
  # Peer repositories whose commits live in another history entirely.
  FOREIGN_REPO_RE <- paste0(
    "stan4bart|bartCause|treatSens|bairrtt|bart-playground|",
    "bart-comp-efficiency|pymc[-_]bart|bartrs|compat branch|dbarts-1\\.0"
  )

  # Collect every (file, line, token) candidate first, then resolve every
  # distinct token with a single batched `git cat-file --batch-check`
  # call - one subprocess total, not one per hash.
  candFile <- character(0)
  candLine <- integer(0)
  candTok <- character(0)

  for (relFile in hashRel) {
    lines <- getContent(relFile)
    if (identical(lines, NA)) {
      next
    }
    # The scope a hash's own citation lives in: one table row, or the
    # whole paragraph a wrapped sentence sits in.
    para <- cumsum(!nzchar(trimws(lines)))
    grpBlob <- vapply(
      split(lines, para),
      function(b) paste(b, collapse = " "),
      character(1L)
    )
    scopeBlob <- ifelse(
      grepl("^\\s*\\|", lines),
      lines,
      unname(grpBlob[as.character(para)])
    )
    inFence <- FALSE
    for (i in seq_along(lines)) {
      line <- lines[i]
      if (grepl("^\\s*```", line)) {
        inFence <- !inFence
        next
      }
      if (inFence) {
        next
      }
      if (grepl(FOREIGN_REPO_RE, scopeBlob[i], perl = TRUE)) {
        next
      }
      m <- gregexpr(hashRe, line, perl = TRUE)[[1L]]
      if (m[1L] == -1L) {
        next
      }
      lens <- attr(m, "match.length")
      for (k in seq_along(m)) {
        tok <- substr(line, m[k], m[k] + lens[k] - 1L)
        if (!grepl("[a-fA-F]", tok)) {
          next
        }
        ctx <- substr(
          line,
          max(1L, m[k] - 80L),
          min(nchar(line), m[k] + lens[k] + 80L)
        )
        if (grepl("pre-rebase", ctx, ignore.case = TRUE)) {
          next
        }
        # Not a commit reference at all: a hex-looking token that is
        # itself part of a filesystem path (e.g. a ".claude/jobs/<id>/"
        # background-job directory) or directly named as a job id.
        charBefore <- if (m[k] > 1L) substr(line, m[k] - 1L, m[k] - 1L) else ""
        charAfter <- if (m[k] + lens[k] <= nchar(line)) {
          substr(line, m[k] + lens[k], m[k] + lens[k])
        } else {
          ""
        }
        precedingWord <- substr(line, max(1L, m[k] - 12L), m[k] - 1L)
        if (
          identical(charBefore, "/") ||
            identical(charAfter, "/") ||
            grepl("\\bjob\\s+$", precedingWord, ignore.case = TRUE)
        ) {
          next
        }
        candFile <- c(candFile, relFile)
        candLine <- c(candLine, i)
        candTok <- c(candTok, tok)
      }
    }
  }

  nHashChecked <- length(candTok)
  if (nHashChecked > 0L) {
    uniqueToks <- unique(candTok)
    # Plain --batch-check (default format, "<sha40> <type> <size>" or
    # "<object> missing") - a custom --batch-check=FORMAT string with
    # parens defeats system2()'s shell quoting on some platforms.
    batchOut <- suppressWarnings(system2(
      "git",
      c("-C", root, "cat-file", "--batch-check"),
      input = paste0(uniqueToks, "^{commit}"),
      stdout = TRUE,
      stderr = FALSE
    ))
    # Each output line is either "<sha40> commit <size>" (resolved) or
    # "<requested-token>^{commit} missing" (echoes the request verbatim).
    okCount <- length(grep("^[0-9a-fA-F]{40} commit [0-9]+$", batchOut))
    missingToks <- sub(
      "\\^\\{commit\\} missing$",
      "",
      grep("\\^\\{commit\\} missing$", batchOut, value = TRUE)
    )
    resolved <- !(uniqueToks %in% missingToks)
    # A malformed/short batch-check response (should not happen) leaves
    # every token unresolved rather than silently passing.
    if (okCount + length(missingToks) != length(uniqueToks)) {
      resolved[] <- FALSE
    }
    names(resolved) <- uniqueToks
    for (idx in which(!resolved[candTok])) {
      report(sprintf(
        "%s:%d: commit hash '%s' does not resolve",
        candFile[idx],
        candLine[idx],
        candTok[idx]
      ))
    }
  }
}

# ---------------------------------------------------------------------
# Part 5: INDEX.md status labels vs. each doc's own Status: line
# ---------------------------------------------------------------------
#
# For every docs/design/*.md carrying a literal "Status:" line in its
# first 12 lines, the leading ALL-CAPS status phrase INDEX.md gives that
# doc (e.g. "MIXED" out of "MIXED (research survey; ...)") must be
# confirmable - every token of that phrase found, case-insensitively, in
# the doc's own first-12-lines block. A miss is a real drift: either the
# doc's header was never updated after a verdict changed, or INDEX.md's
# row was. Docs without a Status: line are counted and skipped (INDEX.md
# itself documents some standing-reference docs as carrying none).

nStatusChecked <- 0L

extractStatusTokens <- function(text) {
  text <- trimws(text)
  tokens <- character(0)
  tokRe <- "^[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*"
  connRe <- "^(,\\s+|\\s*-\\s*|\\s*/\\s*|\\s+)"
  repeat {
    m <- regmatches(text, regexpr(tokRe, text, perl = TRUE))
    if (length(m) == 0L || !nzchar(m)) {
      break
    }
    tokens <- c(tokens, m)
    rest <- substring(text, nchar(m) + 1L)
    cm <- regmatches(rest, regexpr(connRe, rest, perl = TRUE))
    if (length(cm) == 0L || !nzchar(cm)) {
      break
    }
    afterConn <- substring(rest, nchar(cm) + 1L)
    nt <- regmatches(afterConn, regexpr(tokRe, afterConn, perl = TRUE))
    if (length(nt) == 0L || !nzchar(nt)) {
      break
    }
    text <- afterConn
  }
  tokens
}

indexFile2 <- p("docs/design/INDEX.md")
if (file.exists(indexFile2)) {
  idxLines <- readLines(indexFile2, warn = FALSE)
  rowRe <- "^\\|\\s*([A-Za-z0-9_.+-]+\\.md)\\s*\\|\\s*([^|]*?)\\s*\\|"
  for (line in idxLines) {
    rm <- regmatches(line, regexec(rowRe, line, perl = TRUE))[[1L]]
    if (length(rm) < 3L) {
      next
    }
    fname <- rm[2L]
    statusCell <- rm[3L]
    if (identical(fname, "INDEX.md")) {
      next
    }
    docPath <- file.path("docs/design", fname)
    docLines <- getContent(docPath)
    if (identical(docLines, NA)) {
      next # already flagged by Part 1 (INDEX completeness)
    }
    headLines <- docLines[seq_len(min(12L, length(docLines)))]
    if (length(grep("^Status:", headLines)) == 0L) {
      next # no Status: line - counted (via nStatusChecked not being
      # incremented) and skipped, per spec
    }
    nStatusChecked <- nStatusChecked + 1L
    tokens <- extractStatusTokens(statusCell)
    if (length(tokens) == 0L) {
      next
    }
    headBlob <- tolower(paste(headLines, collapse = " "))
    missing <- tokens[
      !vapply(
        tokens,
        function(t) grepl(tolower(t), headBlob, fixed = TRUE),
        logical(1L)
      )
    ]
    if (length(missing) > 0L) {
      report(sprintf(
        "INDEX: %s status '%s' (INDEX.md) not confirmed by %s's own Status: line (missing %s)",
        fname,
        paste(tokens, collapse = " "),
        docPath,
        paste(missing, collapse = ", ")
      ))
    }
  }
}

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

if (length(warnFind) > 0L) {
  cat(paste0("WARN: ", warnFind), sep = "\n")
}
if (length(find) > 0L) {
  cat(paste0("FAIL: ", find), sep = "\n")
  cat(sprintf(
    "check-doc-freshness: %d failure(s), %d warning(s)\n",
    length(find),
    length(warnFind)
  ))
  quit(status = 1L, save = "no")
}
cat(sprintf(
  paste0(
    "check-doc-freshness: OK (%d design docs, %d plan docs indexed; ",
    "%d anchors, %d symbols, %d scenario-count claims verified; ",
    "%d at-line anchors [%d strict, %d advisory], %d warning(s), ",
    "%d ambiguous bare filenames skipped; %d commit hashes verified; ",
    "%d status labels checked)\n"
  ),
  nDesign,
  nPlans,
  nAnchors,
  nSymbols,
  nCounts,
  nAtLineAnchors,
  nAtLineStrict,
  nAtLineAdvisory,
  length(warnFind),
  nAmbiguousAnchor,
  nHashChecked,
  nStatusChecked
))
quit(status = 0L, save = "no")
