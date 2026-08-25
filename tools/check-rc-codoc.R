#!/usr/bin/env Rscript

# Guards dbartsSampler's reference-class methods against usage drift.
#
# tools::codoc() walks registered S3/S4 methods and ordinary closures, not
# a setRefClass generator's methods = list(...): a formal renamed on the
# live generator (R/dbarts.R's dbartsSampler <- setRefClass(...)) never
# fails R CMD check even though a caller who follows the documented
# argument name gets R's own "unused argument" error. This script parses
# both sides directly - the generator's methods list in R/dbarts.R and the
# \S4method{NAME}{dbartsSampler}(...) entries in \usage{} of
# man/dbartsSampler-class.Rd - and compares, per documented method, the
# argument names (order-sensitive) and, where both sides give one, the
# default expression.
#
# Only methods that ARE documented (appear as a \S4method{...} usage
# entry) are checked; a method with no usage entry is out of scope here.
#
# Usage: Rscript tools/check-rc-codoc.R [pkg-root]
# Exit status: 0 if every documented method's usage matches its live
# formals, 1 if any mismatch (failures print with a leading "FAIL:").

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[1L] else "."
p <- function(...) file.path(root, ...)

find <- character(0L)
report <- function(msg) find <<- c(find, msg)

deparseOneLine <- function(expr) paste(trimws(deparse(expr)), collapse = " ")

# ---------------------------------------------------------------------
# Part 1: RC methods live in R/dbarts.R's dbartsSampler <- setRefClass(...)
# ---------------------------------------------------------------------

rFile <- p("R/dbarts.R")
exprs <- parse(file = rFile, keep.source = FALSE)

generatorCall <- NULL
for (e in exprs) {
  if (
    is.call(e) &&
      identical(e[[1]], as.symbol("<-")) &&
      identical(e[[2]], as.symbol("dbartsSampler"))
  ) {
    generatorCall <- e[[3]]
    break
  }
}
if (is.null(generatorCall)) {
  cat(sprintf(
    "FAIL: could not find 'dbartsSampler <- setRefClass(...)' in %s\n",
    rFile
  ))
  quit(status = 1L, save = "no")
}

methodsIdx <- which(names(generatorCall) == "methods")
if (length(methodsIdx) != 1L) {
  cat("FAIL: setRefClass(...) call has no single 'methods =' argument\n")
  quit(status = 1L, save = "no")
}
methodsCall <- generatorCall[[methodsIdx]]
methodNames <- names(methodsCall)[-1L]

liveFormals <- list()
for (i in seq_along(methodNames)) {
  fn <- methodsCall[[i + 1L]]
  if (!is.call(fn) || !identical(fn[[1L]], as.symbol("function"))) {
    next
  }
  fmls <- fn[[2L]]
  nms <- names(fmls)
  if (is.null(nms)) {
    nms <- character(0L)
  }
  defaults <- vapply(
    as.list(fmls),
    function(x) {
      if (identical(x, quote(expr = ))) NA_character_ else deparseOneLine(x)
    },
    character(1L),
    USE.NAMES = FALSE
  )
  liveFormals[[methodNames[i]]] <- list(names = nms, defaults = defaults)
}

# ---------------------------------------------------------------------
# Part 2: \S4method{NAME}{dbartsSampler}(...) entries in \usage{}
# ---------------------------------------------------------------------

rdFile <- p("man/dbartsSampler-class.Rd")
rdText <- paste(readLines(rdFile, warn = FALSE), collapse = "\n")

usageMatch <- regexpr("\\\\usage\\{", rdText)
if (usageMatch < 0) {
  cat(sprintf("FAIL: no \\usage{ block found in %s\n", rdFile))
  quit(status = 1L, save = "no")
}
chars <- strsplit(rdText, "", fixed = TRUE)[[1L]]
startBrace <- usageMatch + attr(usageMatch, "match.length") - 1L
depth <- 0L
endBrace <- NA_integer_
for (i in seq(startBrace, length(chars))) {
  if (chars[i] == "{") {
    depth <- depth + 1L
  } else if (chars[i] == "}") {
    depth <- depth - 1L
    if (depth == 0L) {
      endBrace <- i
      break
    }
  }
}
if (is.na(endBrace)) {
  cat(sprintf("FAIL: unbalanced \\usage{ block in %s\n", rdFile))
  quit(status = 1L, save = "no")
}
usageText <- substr(rdText, startBrace + 1L, endBrace - 1L)
uchars <- strsplit(usageText, "", fixed = TRUE)[[1L]]

# Splits an argument-list string on its top-level commas, respecting
# parenthesis nesting (a default like 'c(nodeHeight = 12, ...)' must not
# be split on its own internal commas).
splitTopLevel <- function(text) {
  text <- trimws(text)
  if (text == "") {
    return(character(0L))
  }
  cs <- strsplit(text, "", fixed = TRUE)[[1L]]
  parenDepth <- 0L
  pieces <- character(0L)
  cur <- character(0L)
  for (ch in cs) {
    if (ch == "(") {
      parenDepth <- parenDepth + 1L
    }
    if (ch == ")") {
      parenDepth <- parenDepth - 1L
    }
    if (ch == "," && parenDepth == 0L) {
      pieces <- c(pieces, paste(cur, collapse = ""))
      cur <- character(0L)
    } else {
      cur <- c(cur, ch)
    }
  }
  pieces <- c(pieces, paste(cur, collapse = ""))
  trimws(pieces)
}

methodPattern <- "\\\\S4method\\{([A-Za-z0-9_.]+)\\}\\{dbartsSampler\\}\\("
m <- gregexpr(methodPattern, usageText, perl = TRUE)[[1L]]
if (m[1L] == -1L) {
  cat("FAIL: no \\S4method{...}{dbartsSampler}(...) usage entries found\n")
  quit(status = 1L, save = "no")
}
matchLens <- attr(m, "match.length")
capStarts <- attr(m, "capture.start")
capLens <- attr(m, "capture.length")

docFormals <- list()
for (k in seq_along(m)) {
  name <- substr(usageText, capStarts[k], capStarts[k] + capLens[k] - 1L)
  parenStart <- m[k] + matchLens[k] - 1L
  depth <- 0L
  parenEnd <- NA_integer_
  for (i in seq(parenStart, length(uchars))) {
    if (uchars[i] == "(") {
      depth <- depth + 1L
    } else if (uchars[i] == ")") {
      depth <- depth - 1L
      if (depth == 0L) {
        parenEnd <- i
        break
      }
    }
  }
  if (is.na(parenEnd)) {
    report(sprintf("unbalanced parens after \\S4method{%s}", name))
    next
  }
  argText <- gsub(
    "[ \t\n]+",
    " ",
    substr(usageText, parenStart + 1L, parenEnd - 1L)
  )
  pieces <- splitTopLevel(argText)
  argNames <- character(length(pieces))
  argDefaults <- character(length(pieces))
  for (j in seq_along(pieces)) {
    piece <- pieces[j]
    eqPos <- regexpr("=", piece, fixed = TRUE)
    if (piece == "...") {
      argNames[j] <- "..."
      argDefaults[j] <- NA_character_
    } else if (eqPos < 0) {
      argNames[j] <- trimws(piece)
      argDefaults[j] <- NA_character_
    } else {
      argNames[j] <- trimws(substr(piece, 1L, eqPos - 1L))
      argDefaults[j] <- trimws(substr(piece, eqPos + 1L, nchar(piece)))
    }
  }
  docFormals[[name]] <- list(names = argNames, defaults = argDefaults)
}

# ---------------------------------------------------------------------
# Part 3: compare
# ---------------------------------------------------------------------

for (name in names(docFormals)) {
  live <- liveFormals[[name]]
  doc <- docFormals[[name]]
  if (is.null(live)) {
    report(sprintf(
      "\\S4method{%s} is documented but no such method exists on dbartsSampler",
      name
    ))
    next
  }
  if (!identical(doc$names, live$names)) {
    report(sprintf(
      "%s: usage names (%s) != live formals (%s)",
      name,
      paste(doc$names, collapse = ", "),
      paste(live$names, collapse = ", ")
    ))
    next
  }
  for (j in seq_along(doc$names)) {
    dDef <- doc$defaults[j]
    lDef <- live$defaults[j]
    if (
      !is.na(dDef) &&
        !is.na(lDef) &&
        !identical(gsub(" ", "", dDef), gsub(" ", "", lDef))
    ) {
      report(sprintf(
        "%s: usage default for '%s' (%s) != live default (%s)",
        name,
        doc$names[j],
        dDef,
        lDef
      ))
    }
  }
}

if (length(find) > 0L) {
  cat(paste0("FAIL: ", find), sep = "\n")
  cat(sprintf("check-rc-codoc: %d mismatch(es)\n", length(find)))
  quit(status = 1L, save = "no")
}
cat(sprintf(
  "check-rc-codoc: OK (%d documented dbartsSampler methods checked)\n",
  length(docFormals)
))
quit(status = 0L, save = "no")
