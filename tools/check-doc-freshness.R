#!/usr/bin/env Rscript

# Guards the documentation's citations and manifests against silent drift.
#
# THE CITE GRAMMAR. Documentation cites code by SYMBOL, never by line
# number. Every citation is delimited by double brackets so one regex
# finds all of them, and the path is separated from what it names by
# "#" so a C++ "::" qualification never collides with the path:
#
#   [[path#symbol]]        CURRENT STATE. The path resolves (alias,
#                          repo path, or a basename unique in the
#                          tracked inventory); the symbol, split on
#                          "::" and "$", must occur in that file with
#                          every component a whole-word token. A
#                          trailing "()" is stripped. Several names in
#                          one cite are separated by commas.
#   [[path#"fragment"]]    VERBATIM. The quoted text must occur in the
#                          file as a literal substring. This is how a
#                          file with no symbols to name - a tinytest
#                          script - is cited; tests/cpp cites name their
#                          test function instead.
#   [[doc.md#Heading]]     DOC TO DOC. A ".md" target is always a
#                          heading cite: the text must appear in one of
#                          the target's markdown headings.
#   [[path:line@sha]]      HISTORY. "line" or "a-b"; the sha is at least
#   [[path:a-b@sha]]       8 hex digits and must be an ancestor of HEAD.
#                          Ancestry is all that is checked - the file
#                          and the line belong to that commit, not to
#                          the working tree. This is the ONLY form in
#                          which a line number may appear, and it is
#                          what landing notes and a plan's "Landing"
#                          section use.
#   retired: [[...]]       SKIPPED. The cite names something that is
#                          gone; the surrounding prose must say so.
#
# Backticks around a cite are optional and accepted either way. Fenced
# code blocks are skipped whole: they quote commands and output, not
# citations.
#
# RESIDUE SELF-CHECK. Every covered line is re-scanned, with the cites
# masked out, for a path-shaped token followed by ":digits", for a bare
# ":digits" continuation, and for the unbracketed "path (symbol)" form.
# Any hit is an unconverted citation and fails: no cite can slip through
# unparsed by being written in a shape the cite regex does not see. A
# "[[" left unclosed on its own line fails for the same reason - the
# cite regex is per-line, so a cite hard-wrapped across a line break
# would otherwise be invisible to every check here.
#
# COVERAGE. Everything under docs/, the top-level README.md, man/*.Rd
# and vignettes/*.Rmd. The last two carry no cites today and are guarded
# so that drift cannot start there.
#
# The other checks, none of which read a citation:
#
# 1. INDEX COMPLETENESS: every docs/design/*.md file is listed in
#    docs/design/INDEX.md, and every name INDEX.md lists exists on disk
#    (no phantom rows). docs/plans/INDEX.md and docs/plans/*.md get the
#    same check when the plans index exists.
#
# 2. SCENARIO COUNTS: docs/design/feature-matrix.md states how many
#    scenarios each equivalence baseline carries; the count is recomputed
#    from the scenario names the recording harness itself defines
#    (benchmarks/R/equivalence.R, bcf-equivalence.R,
#    multinomial-equivalence.R) and a mismatch is flagged.
#
# 3. HASH RESOLVABILITY: every 7-12 character hex token that looks like a
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
# 4. INDEX STATUS: for every docs/design/*.md carrying a literal
#    "Status:" line in its first 12 lines, the leading status phrase
#    INDEX.md gives that doc must be confirmable, token by token,
#    case-insensitively, somewhere in the doc's own first-12-lines
#    block - catching either side going stale after a verdict changes.
#    Docs without a Status: line are counted and skipped.
#
# This is a DEAD-REFERENCE and DRIFTED-COUNT detector, not a semantic
# reviewer: it never re-adjudicates a SHIPPED/REFUSED/etc. cell VALUE,
# and a symbol cite asks only whether the name still occurs in the named
# file, not whether the sentence around it is still true.
#
# Usage: Rscript tools/check-doc-freshness.R [pkg-root]
#                                            [--only <doc path>]...
# --only restricts the cite and residue scan to the named documents and
# skips the four whole-repo checks above.
# Exit status: 0 if every check passes, 1 if any fails. Failures print
# with a leading "FAIL:" so they are easy to grep, and the summary
# breaks the count down by file.

args <- commandArgs(trailingOnly = TRUE)
onlyDocs <- character(0L)
positional <- character(0L)
i <- 1L
while (i <= length(args)) {
  if (args[i] %in% c("--only", "--file")) {
    if (i == length(args)) {
      stop("--only needs a document path")
    }
    onlyDocs <- c(onlyDocs, args[i + 1L])
    i <- i + 2L
  } else {
    positional <- c(positional, args[i])
    i <- i + 1L
  }
}
root <- if (length(positional) >= 1L) positional[1L] else "."
p <- function(...) file.path(root, ...)
restricted <- length(onlyDocs) > 0L

findFile <- character(0L)
findMsg <- character(0L)
report <- function(file, msg) {
  findFile <<- c(findFile, file)
  findMsg <<- c(findMsg, msg)
}

# ---------------------------------------------------------------------
# Shared: path aliases, repo file inventory, per-file content cache
# ---------------------------------------------------------------------

# Path aliases, mirroring the "Path aliases used in cites" block at the
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
  bart.Rd = "man/bart.Rd",
  bart2.Rd = "man/bart2.Rd",
  bartcoreHandle.R = "inst/common/bartcoreHandle.R"
)

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

fileCache <- new.env(parent = emptyenv())
getContent <- function(relPath) {
  if (!exists(relPath, envir = fileCache, inherits = FALSE)) {
    full <- p(relPath)
    val <- if (file.exists(full)) readLines(full, warn = FALSE) else NA
    assign(relPath, val, envir = fileCache)
  }
  get(relPath, envir = fileCache, inherits = FALSE)
}
fileExistsCached <- function(relPath) !identical(getContent(relPath), NA)

# Resolution outcome: the repo-relative path, or one of the two refusal
# codes, which the caller turns into a failure naming the reason.
UNRESOLVED <- "\001unresolved"
AMBIGUOUS <- "\001ambiguous"

resolvePath <- function(tok) {
  if (tok %in% names(ALIAS_PATH)) {
    return(unname(ALIAS_PATH[tok]))
  }
  if (tok %in% names(RD_ALIAS)) {
    return(unname(RD_ALIAS[tok]))
  }
  if (tok %in% R_ALIAS_FILES) {
    return(file.path("R", tok))
  }
  if (grepl("/", tok, fixed = TRUE)) {
    return(if (fileExistsCached(tok)) tok else UNRESOLVED)
  }
  hit <- basenameIndex[[tok]]
  if (is.null(hit)) {
    return(if (fileExistsCached(tok)) tok else UNRESOLVED)
  }
  hit <- hit[vapply(hit, fileExistsCached, logical(1L))]
  if (length(hit) == 0L) {
    return(UNRESOLVED)
  }
  if (length(hit) > 1L) {
    return(AMBIGUOUS)
  }
  hit[1L]
}

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
    report(
      indexRel,
      sprintf("INDEX: %s/%s exists on disk but is not listed", dir, f)
    )
  }
  for (f in phantom) {
    report(
      indexRel,
      sprintf("INDEX: %s/%s is listed but does not exist on disk", dir, f)
    )
  }
  length(onDisk)
}

nDesign <- 0L
nPlans <- 0L
if (!restricted) {
  if (file.exists(p("docs/design/INDEX.md"))) {
    nDesign <- checkIndex("docs/design", "docs/design/INDEX.md")
  } else {
    report("docs/design/INDEX.md", "INDEX: docs/design/INDEX.md is missing")
  }
  if (file.exists(p("docs/plans/INDEX.md"))) {
    nPlans <- checkIndex("docs/plans", "docs/plans/INDEX.md")
  }
}

# ---------------------------------------------------------------------
# Part 2: the cite check
# ---------------------------------------------------------------------

coveredDocs <- function() {
  if (restricted) {
    return(onlyDocs)
  }
  docs <- character(0L)
  if (dir.exists(p("docs"))) {
    docs <- file.path(
      "docs",
      list.files(p("docs"), pattern = "\\.md$", recursive = TRUE)
    )
  }
  if (file.exists(p("README.md"))) {
    docs <- c(docs, "README.md")
  }
  if (dir.exists(p("man"))) {
    docs <- c(docs, file.path("man", list.files(p("man"), pattern = "\\.Rd$")))
  }
  if (dir.exists(p("vignettes"))) {
    docs <- c(
      docs,
      file.path("vignettes", list.files(p("vignettes"), pattern = "\\.Rmd$"))
    )
  }
  sort(docs)
}

docsRel <- coveredDocs()

CITE_RE <- paste0(
  "\\[\\[([^\\]#:@]+)",
  "(?:#([^\\]]+)|:([0-9]+)(?:-([0-9]+))?@([0-9a-fA-F]{8,}))",
  "\\]\\]"
)

# Whole-word occurrence. A name carrying a "." is an R-style compound and
# must not match inside a longer dotted name, so its boundaries exclude
# "." as well; a plain identifier uses ordinary word boundaries, so that
# citing `numCuts` is satisfied by `store.numCuts` as much as by the
# declaration.
escapeRe <- function(x) gsub("([^A-Za-z0-9_])", "\\\\\\1", x, perl = TRUE)

tokenPresent <- function(name, content) {
  pat <- if (grepl(".", name, fixed = TRUE)) {
    paste0("(?<![A-Za-z0-9_.])", escapeRe(name), "(?![A-Za-z0-9_.])")
  } else {
    paste0("\\b", escapeRe(name), "\\b")
  }
  any(grepl(pat, content, perl = TRUE))
}

headingTexts <- function(content) {
  hits <- grep("^#{1,6} ", content, value = TRUE)
  gsub("`", "", sub("^#{1,6} +", "", hits))
}

nSymbolCites <- 0L
nQuoteCites <- 0L
nDocCites <- 0L
nHistoryCites <- 0L
nRetiredCites <- 0L
historySha <- character(0L)

# Residue scanners. The path shape is deliberately looser than the cite
# grammar's: a token that merely LOOKS like a path followed by a line
# number is an unconverted citation whether or not it would have
# resolved.
RESIDUE_PATH_RE <- paste0(
  "(?<![A-Za-z0-9_.])",
  "(?:[A-Za-z0-9_.+-]+/)*[A-Za-z0-9_+-]+\\.[A-Za-z]{1,4}:[0-9]+"
)
RESIDUE_ALIAS_RE <- "\\b(?:RIB|CAPI|MOD|CH|FAC|COM|MOV|SAM|TODO):[0-9]+"
RESIDUE_BARE_RE <- "(?<=[\\s`(|/,]):[0-9]{2,}"
RESIDUE_PAREN_RE <- paste0(
  "((?:[A-Za-z0-9_.+-]+/)*[A-Za-z0-9_+-]+\\.[A-Za-z]{1,4}) ",
  "\\(([A-Za-z_.][A-Za-z0-9_.$]*(?:::[A-Za-z0-9_]+)*(?:\\(\\))?)\\)"
)
KNOWN_EXT_RE <- "\\.(?:R|Rd|Rmd|md|hpp|cpp|cc|c|h|py|in|ac|yaml|yml|csv|rds)$"

residueHits <- function(line) {
  out <- character(0L)
  # An opener with no closer on the same line: a cite wrapped over a line
  # break, which the per-line cite regex would never see.
  openers <- gregexpr("\\[\\[", line, perl = TRUE)[[1L]]
  if (openers[1L] != -1L) {
    closers <- gregexpr("\\]\\]", line, perl = TRUE)[[1L]]
    lastClose <- if (closers[1L] == -1L) -1L else max(as.integer(closers))
    if (max(as.integer(openers)) > lastClose) {
      out <- c(out, "[[ with no ]] on the same line")
    }
  }
  for (re in c(RESIDUE_PATH_RE, RESIDUE_ALIAS_RE, RESIDUE_BARE_RE)) {
    m <- regmatches(line, gregexpr(re, line, perl = TRUE))[[1L]]
    out <- c(out, m[nzchar(m)])
  }
  m <- regmatches(line, gregexpr(RESIDUE_PAREN_RE, line, perl = TRUE))[[1L]]
  for (s in m[nzchar(m)]) {
    g <- regmatches(s, regexec(RESIDUE_PAREN_RE, s, perl = TRUE))[[1L]]
    sym <- g[3L]
    # A parenthesised aside is not a symbol cite: it must read as an
    # identifier (mixed case, an underscore, a "::"/"$" qualifier or a
    # trailing "()") and must not itself be a filename or a bare
    # all-caps status word.
    looksSymbol <- grepl("[a-z]", sym) &&
      grepl("[A-Z_$]|::|\\(\\)", sym) &&
      !grepl(KNOWN_EXT_RE, sym, perl = TRUE)
    if (looksSymbol && !identical(resolvePath(g[2L]), UNRESOLVED)) {
      out <- c(out, s)
    }
  }
  unique(out)
}

for (relDoc in docsRel) {
  lines <- getContent(relDoc)
  if (identical(lines, NA)) {
    report(relDoc, sprintf("%s: does not exist", relDoc))
    next
  }
  isMarkdown <- grepl("\\.(md|Rmd)$", relDoc)
  inFence <- FALSE
  for (ln in seq_along(lines)) {
    line <- lines[ln]
    if (isMarkdown && grepl("^\\s*```", line)) {
      inFence <- !inFence
      next
    }
    if (inFence) {
      next
    }

    m <- gregexpr(CITE_RE, line, perl = TRUE)[[1L]]
    masked <- line
    if (m[1L] != -1L) {
      starts <- as.integer(m)
      lens <- attr(m, "match.length")
      capStart <- attr(m, "capture.start")
      capLen <- attr(m, "capture.length")
      grab <- function(k, g) {
        if (capLen[k, g] <= 0L) {
          return(NA_character_)
        }
        substr(line, capStart[k, g], capStart[k, g] + capLen[k, g] - 1L)
      }
      for (k in seq_along(starts)) {
        # `retired:` immediately before the cite (one optional backtick
        # and any spaces between) marks it as deliberately dead.
        lead <- sub("[` ]*$", "", substr(line, 1L, starts[k] - 1L))
        if (grepl("retired:$", lead)) {
          nRetiredCites <- nRetiredCites + 1L
          next
        }
        pathTok <- trimws(grab(k, 1L))
        payload <- grab(k, 2L)
        cite <- substr(line, starts[k], starts[k] + lens[k] - 1L)

        if (is.na(payload)) {
          # History cite: ancestry, and nothing else.
          nHistoryCites <- nHistoryCites + 1L
          historySha <- c(historySha, tolower(grab(k, 5L)))
          next
        }

        payload <- trimws(payload)
        relPath <- resolvePath(pathTok)
        if (identical(relPath, UNRESOLVED)) {
          report(
            relDoc,
            sprintf("%s:%d: %s does not resolve to a file", relDoc, ln, cite)
          )
          next
        }
        if (identical(relPath, AMBIGUOUS)) {
          report(
            relDoc,
            sprintf(
              "%s:%d: %s is an ambiguous basename - spell the path",
              relDoc,
              ln,
              cite
            )
          )
          next
        }
        content <- getContent(relPath)

        if (grepl("^\".*\"$", payload)) {
          nQuoteCites <- nQuoteCites + 1L
          frag <- substr(payload, 2L, nchar(payload) - 1L)
          if (!any(grepl(frag, content, fixed = TRUE))) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s is not a substring of %s",
                relDoc,
                ln,
                cite,
                relPath
              )
            )
          }
          next
        }

        if (grepl("\\.md$", relPath)) {
          nDocCites <- nDocCites + 1L
          heads <- headingTexts(content)
          if (!any(grepl(payload, heads, fixed = TRUE))) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s names no heading in %s",
                relDoc,
                ln,
                cite,
                relPath
              )
            )
          }
          next
        }

        nSymbolCites <- nSymbolCites + 1L
        names <- trimws(strsplit(payload, ",", fixed = TRUE)[[1L]])
        names <- names[nzchar(names)]
        for (nm in names) {
          nm <- sub("\\(\\)$", "", nm)
          parts <- unlist(strsplit(nm, "::|\\$"))
          parts <- parts[nzchar(parts)]
          absent <- parts[
            !vapply(
              parts,
              tokenPresent,
              logical(1L),
              content = content
            )
          ]
          if (length(absent) > 0L) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s - `%s` not found in %s as a whole-word token",
                relDoc,
                ln,
                cite,
                paste(absent, collapse = "`, `"),
                relPath
              )
            )
          }
        }
      }
      for (k in seq_along(starts)) {
        substr(masked, starts[k], starts[k] + lens[k] - 1L) <- strrep(
          " ",
          lens[k]
        )
      }
    }

    for (hit in residueHits(masked)) {
      report(
        relDoc,
        sprintf(
          "%s:%d: unconverted citation '%s' - cite by symbol inside [[ ]]",
          relDoc,
          ln,
          hit
        )
      )
    }
  }
}

# History-cite ancestry, batched over the distinct shas. A shallow
# checkout (CI) has no history to walk, so the check is skipped there.
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

nShaChecked <- 0L
if (length(historySha) > 0L && !isShallow) {
  uniqueSha <- unique(historySha)
  nShaChecked <- length(uniqueSha)
  bad <- uniqueSha[
    !vapply(
      uniqueSha,
      function(sha) {
        status <- suppressWarnings(system2(
          "git",
          c("-C", root, "merge-base", "--is-ancestor", sha, "HEAD"),
          stdout = FALSE,
          stderr = FALSE
        ))
        identical(as.integer(status), 0L)
      },
      logical(1L)
    )
  ]
  for (sha in bad) {
    where <- character(0L)
    for (relDoc in docsRel) {
      lines <- getContent(relDoc)
      if (identical(lines, NA)) {
        next
      }
      hitLines <- grep(sha, lines, ignore.case = TRUE)
      if (length(hitLines) > 0L) {
        where <- c(where, sprintf("%s:%d", relDoc, hitLines[1L]))
      }
    }
    target <- if (length(where) > 0L) where[1L] else "docs"
    report(
      sub(":[0-9]+$", "", target),
      sprintf(
        "%s: history cite sha '%s' is not an ancestor of HEAD",
        target,
        sha
      )
    )
  }
}

# ---------------------------------------------------------------------
# Part 3: recomputed scenario-count claims (feature-matrix footnotes)
# ---------------------------------------------------------------------
#
# Each baseline .rds is recorded by a benchmarks/R/*.R harness that builds
# its scenario list as a sequence of `result$<name> <- ...` assignments
# inside one named function; the matrix cites how many scenarios that
# baseline carries, and this recounts the names the harness itself
# defines.

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

nCounts <- 0L
matrixRel <- "docs/design/feature-matrix.md"
if (!restricted) {
  if (!fileExistsCached(matrixRel)) {
    report(matrixRel, sprintf("MATRIX: %s is missing", matrixRel))
  } else {
    blob <- paste(getContent(matrixRel), collapse = "\n")
    baselinePat <- "`([a-zA-Z0-9_.-]+\\.rds)`\\s*\\((\\d+)(?:\\s+scenarios)?\\)"
    matchedStrs <- regmatches(
      blob,
      gregexpr(baselinePat, blob, perl = TRUE)
    )[[1L]]
    for (s in matchedStrs) {
      g <- regmatches(s, regexec(baselinePat, s, perl = TRUE))[[1L]]
      fname <- g[2L]
      claimed <- as.integer(g[3L])
      src <- Find(function(e) startsWith(fname, e$prefix), SCENARIO_SOURCES)
      if (is.null(src)) {
        next
      }
      nCounts <- nCounts + 1L
      baselineRel <- file.path("benchmarks/baselines", fname)
      if (!fileExistsCached(baselineRel)) {
        report(
          matrixRel,
          sprintf(
            "MATRIX: baseline '%s' cited (%d scenarios) does not exist at %s",
            fname,
            claimed,
            baselineRel
          )
        )
        next
      }
      names <- scenarioNamesIn(src$file, src$fn)
      if (identical(names, NA_character_)) {
        report(
          matrixRel,
          sprintf(
            "MATRIX: cannot locate %s() in %s to recount scenarios for '%s'",
            src$fn,
            src$file,
            fname
          )
        )
        next
      }
      if (length(names) != claimed) {
        report(
          matrixRel,
          sprintf(
            "MATRIX: '%s' claims %d scenarios but %s()'s %s has %d",
            fname,
            claimed,
            src$fn,
            src$file,
            length(names)
          )
        )
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
designFilesRel <- if (dir.exists(p("docs/design"))) {
  list.files(p("docs/design"), pattern = "\\.md$")
} else {
  character(0L)
}

if (restricted) {
  nHashChecked <- 0L
} else if (isShallow) {
  cat("hash check skipped (shallow clone)\n")
} else {
  hashRel <- c(
    file.path("docs/design", designFilesRel),
    "TODO",
    "benchmarks/baselines/MANIFEST"
  )
  hashRe <- "\\b[0-9a-fA-F]{7,12}\\b"
  FOREIGN_REPO_RE <- paste0(
    "stan4bart|bartCause|treatSens|bairrtt|bart-playground|",
    "bart-comp-efficiency|pymc[-_]bart|bartrs|compat branch|dbarts-1\\.0"
  )

  candFile <- character(0L)
  candLine <- integer(0L)
  candTok <- character(0L)

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
        # itself part of a filesystem path or directly named as a job id.
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
      report(
        candFile[idx],
        sprintf(
          "%s:%d: commit hash '%s' does not resolve",
          candFile[idx],
          candLine[idx],
          candTok[idx]
        )
      )
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
  tokens <- character(0L)
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

indexRel <- "docs/design/INDEX.md"
if (!restricted && fileExistsCached(indexRel)) {
  idxLines <- getContent(indexRel)
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
      next
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
      report(
        indexRel,
        sprintf(
          "INDEX: %s status '%s' not confirmed by %s's own Status: line (missing %s)",
          fname,
          paste(tokens, collapse = " "),
          docPath,
          paste(missing, collapse = ", ")
        )
      )
    }
  }
}

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

if (length(findMsg) > 0L) {
  cat(paste0("FAIL: ", findMsg), sep = "\n")
  cat("\nfailures by file:\n")
  byFile <- sort(table(findFile), decreasing = TRUE)
  for (nm in names(byFile)) {
    cat(sprintf("  %5d  %s\n", byFile[[nm]], nm))
  }
  cat(sprintf(
    "check-doc-freshness: %d failure(s) in %d file(s)\n",
    length(findMsg),
    length(byFile)
  ))
  quit(status = 1L, save = "no")
}
cat(sprintf(
  paste0(
    "check-doc-freshness: OK (%d design docs, %d plan docs indexed; ",
    "%d documents scanned; %d symbol, %d verbatim, %d doc-to-doc, ",
    "%d history (%d distinct sha), %d retired cites; ",
    "%d scenario-count claims, %d commit hashes, %d status labels)\n"
  ),
  nDesign,
  nPlans,
  length(docsRel),
  nSymbolCites,
  nQuoteCites,
  nDocCites,
  nHistoryCites,
  nShaChecked,
  nRetiredCites,
  nCounts,
  nHashChecked,
  nStatusChecked
))
quit(status = 0L, save = "no")
