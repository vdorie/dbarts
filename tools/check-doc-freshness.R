#!/usr/bin/env Rscript

# Guards the documentation's citations and manifests against silent drift.
#
# THE CITE GRAMMAR. Documentation cites code by SYMBOL, never by line
# number. Every citation is an ordinary markdown link, so it renders and
# resolves on GitHub: the TARGET says which file, written relative to
# the citing document's own directory, and the TEXT says what in that
# file is meant.
#
#   [`Sym`](../../src/f.hpp)  CURRENT STATE. The target resolves to a
#                          file that exists; the symbol, split on "::"
#                          and "$", must occur in it with every
#                          component a whole-word token. A trailing
#                          "()" is stripped. Several symbols of one
#                          file are several links separated by ", ";
#                          that run is one cite, counted and marked as
#                          one.
#   ["fragment"](f.R)      VERBATIM. The quoted text must occur in the
#                          file as a literal substring. This is how a
#                          file with no symbols to name - a tinytest
#                          script - is cited; tests/cpp cites name their
#                          test function instead. A bracket inside the
#                          fragment is backslash-escaped so the link
#                          still parses.
#   [Heading](doc.md#slug) DOC TO DOC. An ".md" target carrying a "#"
#                          fragment is a heading cite: the fragment must
#                          be GitHub's slug - lowercased, everything but
#                          letters, digits, spaces, hyphens and
#                          underscores dropped, spaces to hyphens - of
#                          one of the target's markdown headings, and
#                          the link text must be that heading. A heading
#                          of the citing document itself is "(#slug)",
#                          with no path.
#   [path:a-b](blob URL)   HISTORY. The URL is
#                          https://github.com/vdorie/dbarts/blob/<sha>/
#                          <path>#L<a>-L<b>, "#L<a>" for a single line,
#                          with the full 40-hex sha; the link text
#                          repeats "path:a-b", the claim the URL makes.
#                          The sha must be an ancestor of HEAD, and the
#                          claim is read AT THAT COMMIT, not against the
#                          working tree: the path must exist in that
#                          commit's tree and the cited line must be
#                          within the file's length there. A file
#                          deleted since therefore still cites cleanly,
#                          while a line number that never existed at the
#                          named commit does not. This is the ONLY form
#                          in which a line number may appear, and it is
#                          what landing notes and a plan's "Landing"
#                          section use.
#   retired: <link>        The named CONSTRUCT is gone; its content
#                          check is skipped. The location must still be
#                          real - the target resolves, and a history
#                          cite's sha, path and line are checked as
#                          usual. The surrounding prose must say the
#                          thing is gone.
#   unresolved: <link>     The LOCATION cannot be established: a history
#                          cite in a frozen record whose target could
#                          not be placed at any candidate commit. The
#                          link must still parse and the sha must still
#                          be an ancestor; the path and line are not
#                          checked. History cites only.
#
# A marker binds only when it sits immediately before the link, so
# "un-retired:" and any other word ending in "retired:" does not disarm
# one. A link whose text is its target's own path or basename is a plain
# FILE LINK, not a cite, and is not checked; a link out of the
# repository is left alone. Any OTHER link into a tracked file fails: it
# reads as a cite and is checked by nothing.
#
# TWO STANDING LIMITS, both deliberate. Fenced code blocks are skipped
# whole, for the cite check as well as the residue scan: a fence is a
# code sample or a transcript, not a citation, so a cite parked inside
# one is never checked. And a symbol's presence is a token search over
# the whole file, so a name occurring only in a comment or a string
# literal satisfies its cite; qualify a name the target file defines
# more than once ("Class::method") so the token that answers is the one
# the sentence means.
#
# RESIDUE SELF-CHECK. Every covered line is re-scanned, with the links
# masked out, for a path-shaped token followed by ":digits", for a bare
# ":digits" continuation, and for the unbracketed "path (symbol)" form.
# Any hit is an unconverted citation and fails: no cite can slip through
# unparsed by being written in a shape the link scanner does not see. A
# surviving "[[path#symbol]]" or "[[path:line@sha]]" of the retired
# double-bracket grammar fails for the same reason, as does a "[[" left
# unclosed on its own line - the scan is per-line, so a cite hard-
# wrapped across a line break would otherwise be invisible to every
# check here - and a "[[...]]" span whose first token resolves to a
# tracked file.
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
#    the rest). A shallow checkout cannot resolve any hash, so one is
#    refused outright rather than silently passing this check.
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
# Shared: repo file inventory, per-file content cache
# ---------------------------------------------------------------------

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
    val <- if (file.exists(full) && !dir.exists(full)) {
      readLines(full, warn = FALSE)
    } else {
      NA
    }
    assign(relPath, val, envir = fileCache)
  }
  get(relPath, envir = fileCache, inherits = FALSE)
}
fileExistsCached <- function(relPath) !identical(getContent(relPath), NA)

# Does a bare token name a tracked file? Answered for the residue scan,
# which reads a path out of prose rather than out of a link: the
# repo-relative path, or one of the two refusal codes.
UNRESOLVED <- "\001unresolved"
AMBIGUOUS <- "\001ambiguous"

resolvePath <- function(tok) {
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

# An inline markdown link. The text may carry backslash-escaped brackets
# (a verbatim fragment containing one) but no bare bracket; the target
# holds no whitespace or parenthesis, which every path and every blob URL
# here satisfies.
LINK_RE <- "\\[((?:[^\\[\\]\\\\]|\\\\.)*)\\]\\(([^() \t]*)\\)"
BLOB_PREFIX <- "https://github.com/vdorie/dbarts/blob/"
BLOB_RE <- paste0(
  "^https://github\\.com/vdorie/dbarts/blob/",
  "([0-9a-f]{40})/([^#]+)#L([0-9]+)(?:-L([0-9]+))?$"
)
# The retired double-bracket cite, kept only so the residue scan can fail
# on one that outlived the conversion to links.
OLD_CITE_RE <- paste0(
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

# GitHub's heading anchor: lowercased, every character but a letter,
# digit, space, hyphen or underscore dropped, then spaces to hyphens.
slugOf <- function(heading) {
  s <- tolower(trimws(heading))
  gsub(" ", "-", gsub("[^A-Za-z0-9 _-]", "", s, perl = TRUE), fixed = TRUE)
}

# A bracket inside a link's text is backslash-escaped so the link parses;
# the cite means the character itself.
unescapeText <- function(text) gsub("\\\\(.)", "\\1", text, perl = TRUE)

# A relative target read from the citing document's own directory,
# collapsed lexically rather than through the filesystem so a target that
# climbs past the repo root is refused rather than silently rebased.
lexicalPath <- function(dir, rel) {
  full <- if (nzchar(dir) && !identical(dir, ".")) {
    paste0(dir, "/", rel)
  } else {
    rel
  }
  out <- character(0L)
  for (part in strsplit(full, "/", fixed = TRUE)[[1L]]) {
    if (!nzchar(part) || identical(part, ".")) {
      next
    }
    if (identical(part, "..")) {
      if (length(out) == 0L) {
        return(NA_character_)
      }
      out <- out[-length(out)]
    } else {
      out <- c(out, part)
    }
  }
  paste(out, collapse = "/")
}

# Code spans are opaque to the link scanner: a "[" inside `back ticks`
# opens nothing, exactly as GitHub renders it. A span may run across a
# line break inside one paragraph, so the pairing is computed over the
# whole paragraph and handed back line by line, each span's characters
# replaced (newlines excepted, so the lines still split apart) by a
# filler no link syntax can use.
maskCodeSpans <- function(lines) {
  if (!any(grepl("`", lines, fixed = TRUE))) {
    return(lines)
  }
  blob <- paste(lines, collapse = "\n")
  m <- gregexpr("`+", blob, perl = TRUE)[[1L]]
  starts <- as.integer(m)
  lens <- attr(m, "match.length")
  chars <- strsplit(blob, "", fixed = TRUE)[[1L]]
  used <- rep(FALSE, length(starts))
  i <- 1L
  while (i <= length(starts)) {
    if (!used[i]) {
      j <- i + 1L
      while (j <= length(starts) && (used[j] || lens[j] != lens[i])) {
        j <- j + 1L
      }
      if (j <= length(starts)) {
        idx <- starts[i]:(starts[j] + lens[j] - 1L)
        chars[idx[chars[idx] != "\n"]] <- "\001"
        used[i] <- TRUE
        used[j] <- TRUE
        i <- j
      }
    }
    i <- i + 1L
  }
  strsplit(paste(chars, collapse = ""), "\n", fixed = TRUE)[[1L]]
}

maskedLinesOf <- function(lines, skip) {
  out <- lines
  n <- length(lines)
  i <- 1L
  while (i <= n) {
    if (skip[i] || !nzchar(trimws(lines[i]))) {
      i <- i + 1L
      next
    }
    j <- i
    while (j < n && !skip[j + 1L] && nzchar(trimws(lines[j + 1L]))) {
      j <- j + 1L
    }
    out[i:j] <- maskCodeSpans(lines[i:j])
    i <- j + 1L
  }
  out
}

# What a link is, before any content is read. "external" and "file" are
# the two non-cites: a link out of the repository, and a link whose text
# is the target file's own name.
classifyLink <- function(relDoc, docDir, text, target) {
  if (grepl("^[a-zA-Z][a-zA-Z0-9+.-]*:", target, perl = TRUE)) {
    g <- regmatches(target, regexec(BLOB_RE, target, perl = TRUE))[[1L]]
    if (length(g) == 0L) {
      if (startsWith(target, BLOB_PREFIX)) {
        return(list(
          kind = "bad",
          msg = "is not a history cite - the URL needs a 40-hex sha and an #Lnn line anchor"
        ))
      }
      return(list(kind = "external"))
    }
    return(list(
      kind = "history",
      sha = g[2L],
      path = g[3L],
      lo = as.integer(g[4L]),
      hi = as.integer(if (nzchar(g[5L])) g[5L] else g[4L])
    ))
  }
  parts <- strsplit(target, "#", fixed = TRUE)[[1L]]
  pathPart <- if (length(parts) == 0L) "" else parts[1L]
  frag <- if (length(parts) > 1L) paste(parts[-1L], collapse = "#") else ""
  relPath <- if (!nzchar(pathPart)) relDoc else lexicalPath(docDir, pathPart)
  if (is.na(relPath) || !nzchar(relPath) || !fileExistsCached(relPath)) {
    return(list(kind = "bad", msg = "does not resolve to a file"))
  }
  if (grepl("\\.md$", relPath) && nzchar(frag)) {
    return(list(kind = "doc", path = relPath, frag = frag))
  }
  bare <- trimws(gsub("`", "", text, fixed = TRUE))
  if (grepl("^`.*`$", text)) {
    if (identical(bare, pathPart) || identical(bare, relPath)) {
      return(list(kind = "file"))
    }
    return(list(kind = "symbol", path = relPath))
  }
  if (grepl("^\".*\"$", text)) {
    return(list(kind = "verbatim", path = relPath))
  }
  if (
    identical(bare, pathPart) ||
      identical(bare, relPath) ||
      identical(bare, basename(relPath))
  ) {
    return(list(kind = "file"))
  }
  list(
    kind = "bad",
    msg = "names a tracked file but is no symbol, fragment, heading or path:line cite"
  )
}

nSymbolCites <- 0L
nQuoteCites <- 0L
nDocCites <- 0L
nHistoryCites <- 0L
nRetiredCites <- 0L
nUnresolvedCites <- 0L
# One record per history cite; verified against the named commits after
# the whole corpus is scanned, so the git calls batch.
histRec <- list()

# Residue scanners. The path shape is deliberately looser than the cite
# grammar's: a token that merely LOOKS like a path followed by a line
# number is an unconverted citation whether or not it would have
# resolved.
RESIDUE_PATH_RE <- paste0(
  "(?<![A-Za-z0-9_.])",
  "(?:[A-Za-z0-9_.+-]+/)*[A-Za-z0-9_+-]+\\.[A-Za-z]{1,4}:[0-9]+"
)
RESIDUE_ALIAS_RE <- "\\bTODO:[0-9]+"
# A wrapped citation can leave a bare ":NNN" as the very first thing on a
# line, with no delimiter before it for the lookbehind to see. Match that
# case too, but not a timestamp, which never starts a line at the colon
# itself. The digit and range groups are possessive so a wrapped range
# cannot satisfy the negative lookahead by backtracking off its end.
RESIDUE_BARE_RE <- paste0(
  "(?:^|(?<=[\\s`(|/,])):[0-9]{2,}+(?:-[0-9]++)?+",
  "(?!@[0-9a-fA-F]{8,})"
)
RESIDUE_PAREN_RE <- paste0(
  "((?:[A-Za-z0-9_.+-]+/)*[A-Za-z0-9_+-]+\\.[A-Za-z]{1,4}) ",
  "\\(([A-Za-z_.][A-Za-z0-9_.$]*(?:::[A-Za-z0-9_]+)*(?:\\(\\))?)\\)"
)
KNOWN_EXT_RE <- "\\.(?:R|Rd|Rmd|md|hpp|cpp|cc|c|h|py|in|ac|yaml|yml|csv|rds)$"

# Bracketed spans left after the links are masked out: flagged only when
# the first token resolves to a tracked file, which no R or C++ bracket
# idiom does.
payloadlessSpans <- function(line) {
  if (!grepl("[[", line, fixed = TRUE)) {
    return(character(0L))
  }
  m <- regmatches(line, gregexpr("\\[\\[[^][]*\\]\\]", line, perl = TRUE))[[1L]]
  m <- m[nzchar(m)]
  if (length(m) == 0L) {
    return(character(0L))
  }
  m <- m[!grepl(OLD_CITE_RE, m, perl = TRUE)]
  if (length(m) == 0L) {
    return(character(0L))
  }
  inner <- substr(m, 3L, nchar(m) - 2L)
  tok <- trimws(sub("[[:space:]#:].*$", "", inner))
  keep <- vapply(
    tok,
    function(t) {
      nzchar(t) && !(resolvePath(t) %in% c(UNRESOLVED, AMBIGUOUS))
    },
    logical(1L)
  )
  unique(m[keep])
}

residueHits <- function(line) {
  out <- character(0L)
  # A double-bracket cite that outlived the conversion to links, and an
  # opener with no closer on the same line: a cite wrapped over a line
  # break, which the per-line scan would never see.
  old <- regmatches(line, gregexpr(OLD_CITE_RE, line, perl = TRUE))[[1L]]
  out <- c(out, old[nzchar(old)])
  openers <- if (grepl("[[", line, fixed = TRUE)) {
    gregexpr("\\[\\[", line, perl = TRUE)[[1L]]
  } else {
    -1L
  }
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
  docDir <- dirname(relDoc)
  skip <- logical(length(lines))
  inFence <- FALSE
  for (ln in seq_along(lines)) {
    if (isMarkdown && grepl("^\\s*```", lines[ln])) {
      skip[ln] <- TRUE
      inFence <- !inFence
      next
    }
    skip[ln] <- inFence
  }
  maskedLines <- maskedLinesOf(lines, skip)

  for (ln in seq_along(lines)) {
    if (skip[ln]) {
      next
    }
    line <- lines[ln]
    m <- gregexpr(LINK_RE, maskedLines[ln], perl = TRUE)[[1L]]
    masked <- line
    if (m[1L] != -1L) {
      starts <- as.integer(m)
      lens <- attr(m, "match.length")
      capStart <- attr(m, "capture.start")
      capLen <- attr(m, "capture.length")
      grab <- function(k, g) {
        if (capLen[k, g] <= 0L) {
          return("")
        }
        substr(line, capStart[k, g], capStart[k, g] + capLen[k, g] - 1L)
      }
      prevEnd <- -1L
      prevTarget <- ""
      prevSymbol <- FALSE
      prevRetired <- FALSE
      for (k in seq_along(starts)) {
        text <- grab(k, 1L)
        target <- grab(k, 2L)
        cite <- substr(line, starts[k], starts[k] + lens[k] - 1L)
        # A marker binds only immediately before the link (one optional
        # backtick and any spaces between) and only as a whole word, so
        # a hyphenated "un-retired:" disarms nothing.
        lead <- sub("[` ]*$", "", substr(line, 1L, starts[k] - 1L))
        markerRe <- "(?:^|[^A-Za-z0-9_-])%s:$"
        isRetired <- grepl(sprintf(markerRe, "retired"), lead, perl = TRUE)
        isUnresolved <- grepl(
          sprintf(markerRe, "unresolved"),
          lead,
          perl = TRUE
        )
        cl <- classifyLink(relDoc, docDir, text, target)

        # Several symbols of one file are written as several links
        # separated by ", ": the run is one cite, and a marker in front
        # of it covers the whole of it.
        inRun <- identical(cl$kind, "symbol") &&
          prevSymbol &&
          identical(target, prevTarget) &&
          grepl(
            "^, *$",
            substr(line, prevEnd + 1L, starts[k] - 1L)
          )
        if (inRun) {
          isRetired <- prevRetired
        }
        prevEnd <- starts[k] + lens[k] - 1L
        prevTarget <- target
        prevSymbol <- identical(cl$kind, "symbol")
        prevRetired <- isRetired

        if (identical(cl$kind, "external")) {
          next
        }
        if (identical(cl$kind, "bad")) {
          report(
            relDoc,
            sprintf("%s:%d: %s %s", relDoc, ln, cite, cl$msg)
          )
          next
        }
        if (identical(cl$kind, "file")) {
          next
        }

        if (identical(cl$kind, "history")) {
          nHistoryCites <- nHistoryCites + 1L
          if (isRetired) {
            nRetiredCites <- nRetiredCites + 1L
          }
          if (isUnresolved) {
            nUnresolvedCites <- nUnresolvedCites + 1L
          }
          claim <- if (cl$lo == cl$hi) {
            sprintf("%s:%d", cl$path, cl$lo)
          } else {
            sprintf("%s:%d-%d", cl$path, cl$lo, cl$hi)
          }
          if (!identical(text, claim)) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s - the link text must read `%s`, what the URL claims",
                relDoc,
                ln,
                cite,
                claim
              )
            )
            next
          }
          histRec[[length(histRec) + 1L]] <- list(
            doc = relDoc,
            line = ln,
            cite = cite,
            path = cl$path,
            lo = cl$lo,
            hi = cl$hi,
            sha = cl$sha,
            skipPath = isUnresolved
          )
          next
        }

        if (isUnresolved) {
          report(
            relDoc,
            sprintf(
              "%s:%d: %s - `unresolved:` marks a history cite only",
              relDoc,
              ln,
              cite
            )
          )
          next
        }
        if (isRetired) {
          # The construct is gone, so its content is not checked; the
          # location it names still has to be real, which classifyLink
          # has already established.
          if (!inRun) {
            nRetiredCites <- nRetiredCites + 1L
          }
          next
        }
        content <- getContent(cl$path)

        if (identical(cl$kind, "verbatim")) {
          nQuoteCites <- nQuoteCites + 1L
          frag <- unescapeText(substr(text, 2L, nchar(text) - 1L))
          if (!any(grepl(frag, content, fixed = TRUE))) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s is not a substring of %s",
                relDoc,
                ln,
                cite,
                cl$path
              )
            )
          }
          next
        }

        if (identical(cl$kind, "doc")) {
          nDocCites <- nDocCites + 1L
          heads <- headingTexts(content)
          named <- heads[slugOf(heads) == cl$frag]
          if (length(named) == 0L) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s names no heading in %s",
                relDoc,
                ln,
                cite,
                cl$path
              )
            )
          } else if (!(trimws(unescapeText(text)) %in% trimws(named))) {
            report(
              relDoc,
              sprintf(
                "%s:%d: %s - the link text must be the heading it points at",
                relDoc,
                ln,
                cite
              )
            )
          }
          next
        }

        if (!inRun) {
          nSymbolCites <- nSymbolCites + 1L
        }
        names <- trimws(strsplit(text, ",", fixed = TRUE)[[1L]])
        names <- gsub("^`|`$", "", names)
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
                cl$path
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
          "%s:%d: unconverted citation '%s' - cite by symbol in a markdown link",
          relDoc,
          ln,
          hit
        )
      )
    }
    # A bracketed span no link matched, whose first token nonetheless
    # names a tracked file: it reads as a cite and is checked by
    # nothing. R and C++ syntax ("[[1L]]", "[[nodiscard]]") never names
    # one.
    for (span in payloadlessSpans(masked)) {
      report(
        relDoc,
        sprintf(
          "%s:%d: %s names a file but no symbol, fragment or path:line",
          relDoc,
          ln,
          span
        )
      )
    }
  }
}

# ---------------------------------------------------------------------
# History-cite verification
# ---------------------------------------------------------------------
#
# Three questions per cite, batched over the distinct shas and the
# distinct (sha, path) pairs so the whole corpus costs two git calls
# plus one per distinct sha:
#
#   a. is the sha an ancestor of HEAD?
#   b. does the path exist in that commit's tree?
#   c. is the cited line within that blob's length there?
#
# (b) and (c) are what make the history form a claim rather than an
# assertion: without them a line number that never existed at the named
# commit passes, and a converter that guessed the wrong file for a bare
# ":NNN" leaves no trace. An `unresolved:` cite answers (a) only.
#
# The path is the URL's own, read as written and never checked against
# the current tree, since a history cite may legitimately name a path
# that has since been deleted; git at that commit is the judge.

haveGit <- tryCatch(
  identical(
    as.integer(suppressWarnings(system2(
      "git",
      c("-C", root, "rev-parse", "--git-dir"),
      stdout = FALSE,
      stderr = FALSE
    ))),
    0L
  ),
  error = function(e) FALSE
)
isShallow <- haveGit &&
  identical(
    trimws(suppressWarnings(system2(
      "git",
      c("-C", root, "rev-parse", "--is-shallow-repository"),
      stdout = TRUE,
      stderr = FALSE
    ))),
    "true"
  )
historyUsable <- haveGit && !isShallow
if (!historyUsable) {
  report(
    "tools/check-doc-freshness.R",
    paste0(
      if (haveGit) {
        "HISTORY: the checkout is shallow"
      } else {
        "HISTORY: this is not a git checkout"
      },
      " - commit ancestry, the path-at-commit check and the ",
      "commit-hash check all need full history. Check out with ",
      "fetch-depth: 0 (actions/checkout) or run `git fetch --unshallow`."
    )
  )
}

nShaChecked <- 0L
nHistoryPairs <- 0L

# The line count of every distinct (sha, path), from one `git cat-file
# --batch` stream read as bytes: the header names the blob size, so each
# body is sliced out exactly and its newlines counted. A body whose last
# byte is not a newline still ends a line, matching readLines().
blobLineCounts <- function(specs) {
  outFile <- tempfile("blobs")
  on.exit(unlink(outFile), add = TRUE)
  suppressWarnings(system2(
    "git",
    c("-C", root, "cat-file", "--batch"),
    input = specs,
    stdout = outFile,
    stderr = FALSE
  ))
  found <- rep(FALSE, length(specs))
  lineCount <- rep(NA_integer_, length(specs))
  names(found) <- specs
  names(lineCount) <- specs
  if (!file.exists(outFile) || file.size(outFile) == 0) {
    return(list(found = found, lineCount = lineCount))
  }
  bytes <- readBin(outFile, "raw", file.size(outFile))
  NL <- as.raw(10L)
  pos <- 1L
  for (i in seq_along(specs)) {
    if (pos > length(bytes)) {
      break
    }
    window <- bytes[pos:min(pos + 512L, length(bytes))]
    nl <- which(window == NL)[1L]
    if (is.na(nl)) {
      break
    }
    header <- rawToChar(window[seq_len(nl - 1L)])
    pos <- pos + nl
    if (grepl(" (missing|ambiguous)$", header)) {
      next
    }
    size <- suppressWarnings(as.integer(sub("^[0-9a-f]+ [a-z]+ ", "", header)))
    if (is.na(size)) {
      break
    }
    found[i] <- TRUE
    if (size == 0L) {
      lineCount[i] <- 0L
    } else {
      body <- bytes[pos:(pos + size - 1L)]
      n <- sum(body == NL)
      if (body[length(body)] != NL) {
        n <- n + 1L
      }
      lineCount[i] <- n
    }
    pos <- pos + size + 1L
  }
  list(found = found, lineCount = lineCount)
}

if (length(histRec) > 0L && historyUsable) {
  fieldOf <- function(name, mode) {
    vapply(histRec, function(r) r[[name]], vector(mode, 1L))
  }
  hDoc <- fieldOf("doc", "character")
  hLine <- fieldOf("line", "integer")
  hCite <- fieldOf("cite", "character")
  hTok <- fieldOf("path", "character")
  hLo <- fieldOf("lo", "integer")
  hHi <- fieldOf("hi", "integer")
  hSha <- fieldOf("sha", "character")
  hSkip <- fieldOf("skipPath", "logical")

  # --- (a) ancestry, over the distinct shas ---
  # Two git calls, not one per sha: expand every abbreviation to its full
  # oid, then test membership in the set of commits reachable from HEAD,
  # which is exactly "is an ancestor of HEAD".
  uniqueSha <- unique(hSha)
  nShaChecked <- length(uniqueSha)
  expanded <- suppressWarnings(system2(
    "git",
    c("-C", root, "cat-file", "--batch-check"),
    input = paste0(uniqueSha, "^{commit}"),
    stdout = TRUE,
    stderr = FALSE
  ))
  fullOid <- ifelse(
    grepl("^[0-9a-f]{40} commit [0-9]+$", expanded),
    sub(" .*$", "", expanded),
    NA_character_
  )
  reachable <- suppressWarnings(system2(
    "git",
    c("-C", root, "rev-list", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ))
  ancestor <- !is.na(fullOid) & fullOid %in% reachable
  if (length(ancestor) != length(uniqueSha)) {
    ancestor <- rep(FALSE, length(uniqueSha))
  }
  names(ancestor) <- uniqueSha
  for (i in which(!ancestor[hSha])) {
    report(
      hDoc[i],
      sprintf(
        "%s:%d: %s - sha '%s' is not an ancestor of HEAD",
        hDoc[i],
        hLine[i],
        hCite[i],
        hSha[i]
      )
    )
  }

  # --- (b) and (c), over the distinct (sha, path) pairs ---
  hPath <- hTok
  checkable <- unname(ancestor[hSha]) & !hSkip
  if (any(checkable)) {
    specs <- unique(paste0(hSha[checkable], ":", hPath[checkable]))
    nHistoryPairs <- length(specs)
    blobs <- blobLineCounts(specs)
    key <- paste0(hSha, ":", hPath)
    for (i in which(checkable)) {
      k <- key[i]
      if (!blobs$found[[k]]) {
        report(
          hDoc[i],
          sprintf(
            "%s:%d: %s - %s does not exist at %s",
            hDoc[i],
            hLine[i],
            hCite[i],
            hPath[i],
            hSha[i]
          )
        )
      } else if (hHi[i] > blobs$lineCount[[k]]) {
        report(
          hDoc[i],
          sprintf(
            "%s:%d: %s - %s has %d lines at %s",
            hDoc[i],
            hLine[i],
            hCite[i],
            hPath[i],
            blobs$lineCount[[k]],
            hSha[i]
          )
        )
      }
    }
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
# A shallow checkout cannot resolve any hash; that is refused above,
# with the rest of the history checks, rather than skipped here.

nHashChecked <- 0L
designFilesRel <- if (dir.exists(p("docs/design"))) {
  list.files(p("docs/design"), pattern = "\\.md$")
} else {
  character(0L)
}

if (restricted || !historyUsable) {
  nHashChecked <- 0L
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
    "%d history (%d distinct sha, %d distinct sha:path), ",
    "%d retired, %d unresolved cites; ",
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
  nHistoryPairs,
  nRetiredCites,
  nUnresolvedCites,
  nCounts,
  nHashChecked,
  nStatusChecked
))
quit(status = 0L, save = "no")
