#!/usr/bin/env Rscript

# Executable composition matrix: an oracle for docs/design/feature-matrix.md.
# Every claim the matrix makes is turned into one small composition that is
# actually attempted, and the outcome is classified as CONSTRUCTS (runs without
# error), REFUSES (a named stop()), ERRORS (a generic R/engine failure nothing
# chose to name) or HARNESS (the runner's own fixture or argument was wrong).
# The outcome is then diffed against what the matrix says:
#
#   S  must CONSTRUCT, and where a read-back is cheap the effect must show:
#      latents non-NULL and one per row, sigma read back by getSigmas(), the
#      active-row mask on the sampler, installed trees equal to the donor's,
#      weights, response, offset, predictors and test predictors read back
#      from the data object, a named leaf sd changing the leaf prior
#   R  must REFUSE, and the message must match the package's own refusal text
#      pinned here for the cell; any other error, or a REFUSES with other
#      text (a typo in a fixture, say), disagrees
#   M  must not CONSTRUCT, with the same text rule; a bare R "unused argument"
#      is accepted only where the cell is unreachable for want of a formal
#   -  unchecked, except: unit-scale (footnote f7: updateScale = TRUE is
#      ignored, so it must construct) and latents (gaussian and hetero
#      gaussian report NULL)
#
# Two kinds of claim are read. The two tables (section 1, five capability
# columns; section 2, the xbart() and flat C columns) are parsed out of the
# file cell by cell, qualifiers included, and checked against pins: the exact
# row and column sets, the number of non-"-" alternatives per table and the
# number of checks per part. A cell, qualifier or count this runner has no
# fixture or pin for stops the run, so a cell cannot leave the matrix, or a
# row join it, unnoticed; change a pin deliberately when the matrix changes.
# The prose rules (construction reach, the universal capabilities, zero
# weights, variance forest, DART, warm start and grow-from-root, whole-data
# setData) have no table, so each is a hand-written probe anchored to a phrase
# that must still be in the matrix.
#
# Needs: the installed dbarts (the flat C consumer is compiled from its
# inst/tinytest/capi sources, so a C compiler and make), survival, and the
# matrix file two directories up from this script. The run stops, and says so,
# when the consumer cannot be built.
#
# Does NOT check: statistical correctness of any draw; that updateScale or
# growFromRoot changed anything (they are checked to run); that a variance
# forest's own trees are not DART (only that a DART mean forest composes);
# per-family flat C buffers beyond the dimensions and one plain run; any cell
# the matrix marks "-" other than the two kinds above; the Gaps table.
#
# Lives in benchmarks/R/ rather than tools/: it is a local correctness gate
# against a design document, the same role equivalence.R and bench-sampler.R
# play against a baseline, not a build-freshness check. It runs in no
# workflow; run it by hand at the release candidate.
#
# Every fixture is tiny (n = 40, a handful of trees, 0 burn-in, 1 sample,
# n.threads = 1L throughout) and reseeded per check, so a run is deterministic
# and fast.
#
# Usage: Rscript benchmarks/R/composition-matrix.R [--verbose] [--flip=<text>]
# Exit status is nonzero iff at least one disagreement or harness failure was
# recorded. --verbose prints every check. --flip=<text> inverts the
# expectation of every check whose label contains <text>, to show that a
# check can fail.

suppressPackageStartupMessages(library(dbarts))

cmdArgs <- commandArgs()
scriptDir <- dirname(sub(
  "--file=",
  "",
  grep("--file=", cmdArgs, value = TRUE)
))
flipText <- sub("^--flip=", "", grep("^--flip=", cmdArgs, value = TRUE))
flipText <- if (length(flipText) == 1L) flipText else ""
verbose <- "--verbose" %in% cmdArgs
matrixPath <- file.path(
  scriptDir,
  "..",
  "..",
  "docs",
  "design",
  "feature-matrix.md"
)

## ---- 0. Pins -----------------------------------------------------------------

# Deliberate: update these when the matrix, or the set of checks, changes on
# purpose. nonDash counts the non-"-" alternatives parsed per table.
pinNonDash <- c(table1 = 51L, table2 = 24L)
pinChecks <- c(table1 = 61L, table2 = 27L, prose = 181L)

pinStop <- function(what, got, pinned) {
  stop(
    "composition-matrix: ",
    what,
    " is ",
    got,
    ", pinned at ",
    pinned,
    " - if the matrix or the check set changed on purpose, update the pin ",
    "deliberately; otherwise a cell has left (or joined) the matrix"
  )
}

## ---- 1. Parse the matrix's two tables --------------------------------------

rowKeys <- c(
  "gaussian",
  "student",
  "probit",
  "logistic",
  "ordinal",
  "nbinom",
  "multinom",
  "aft",
  "hazard",
  "hurdle",
  "bcf",
  "hetero"
)
# transcribed once from the doc's table headers; the "start" attribute finds
# the table, the names are dispatch keys of its columns - the cell VALUES are
# read out of the file
tables <- list(
  table1 = structure(
    c(
      setWeights = "case weights",
      setSigma = "sigma",
      getLatents = "latents",
      updateScale = "unit-scale transform",
      testSurface = "test fits defined"
    ),
    start = "| model | case weights |"
  ),
  table2 = structure(
    c(xbart = "`xbart()`", flatc = "flat C `dbarts.h`"),
    start = "| model | `xbart()` |"
  )
)

# a cell is one or more alternatives split by " / "; each is a token (S, R, M,
# -) with an optional qualifier after a comma, or the whole cell "as <row>",
# which copies that row's cell
parseCell <- function(text) {
  alts <- strsplit(text, " / ", fixed = TRUE)[[1L]]
  lapply(alts, function(a) {
    a <- gsub("\\[[^]]*\\]\\([^)]*\\)", "", a) # cites
    a <- gsub("\\[f[0-9]+\\]", "", a) # footnote marks
    a <- gsub("\\s+", " ", trimws(a))
    a <- sub(" ,", ",", a, fixed = TRUE)
    a <- sub("[ ,]+$", "", a)
    if (grepl("^as [a-z]+$", a)) {
      return(list(token = "as", qualifier = sub("^as ", "", a)))
    }
    m <- regmatches(a, regexec("^(S|R|M|-)(?:, (.+))?$", a))[[1L]]
    if (length(m) == 0L) {
      stop(
        "composition-matrix: cell '",
        text,
        "' has no legend token - matrix vocabulary needs a human look"
      )
    }
    list(token = m[2L], qualifier = if (nzchar(m[3L])) m[3L] else "")
  })
}

parseTables <- function(lines) {
  claims <- list()
  for (tab in names(tables)) {
    cols <- tables[[tab]]
    start <- attr(cols, "start")
    hdr <- which(startsWith(trimws(lines), start))
    if (length(hdr) != 1L) {
      stop("composition-matrix: no unique table header '", start, "'")
    }
    headerCells <- trimws(strsplit(lines[hdr], "|", fixed = TRUE)[[1L]][-1L])
    if (!identical(headerCells[-1L], as.character(cols))) {
      stop(
        "composition-matrix: columns under '",
        start,
        "' read '",
        paste(headerCells[-1L], collapse = "' '"),
        "', expected '",
        paste(cols, collapse = "' '"),
        "' - column drift"
      )
    }
    # the rows: every line of the table, in both directions
    last <- hdr + 1L
    while (
      last + 1L <= length(lines) && startsWith(trimws(lines[last + 1L]), "|")
    ) {
      last <- last + 1L
    }
    body <- lines[seq(hdr + 2L, length.out = last - hdr - 1L)]
    keys <- trimws(vapply(
      strsplit(body, "|", fixed = TRUE),
      function(p) p[2L],
      ""
    ))
    if (!identical(keys, rowKeys)) {
      stop(
        "composition-matrix: rows under '",
        start,
        "' read '",
        paste(keys, collapse = "' '"),
        "', expected '",
        paste(rowKeys, collapse = "' '"),
        "' - row drift"
      )
    }
    parsed <- list()
    for (i in seq_along(rowKeys)) {
      raw <- trimws(strsplit(body[i], "|", fixed = TRUE)[[1L]][-1L])
      if (length(raw) != length(cols) + 1L) {
        stop(
          "composition-matrix: row '",
          rowKeys[i],
          "' under '",
          start,
          "' has ",
          length(raw),
          " cells, expected ",
          length(cols) + 1L
        )
      }
      parsed[[rowKeys[i]]] <- lapply(raw[-1L], parseCell)
    }
    n <- 0L
    for (family in rowKeys) {
      for (j in seq_along(cols)) {
        alts <- parsed[[family]][[j]]
        if (alts[[1L]]$token == "as") {
          alts <- parsed[[alts[[1L]]$qualifier]][[j]]
        }
        for (alt in alts) {
          claims[[length(claims) + 1L]] <- list(
            table = tab,
            family = family,
            capability = names(cols)[j],
            token = alt$token,
            qualifier = alt$qualifier
          )
          n <- n + (alt$token != "-")
        }
      }
    }
    if (n != pinNonDash[[tab]]) {
      pinStop(paste0("the non-'-' cell count of ", tab), n, pinNonDash[[tab]])
    }
  }
  claims
}

## ---- 2. Shared fixture data -----------------------------------------------

n <- 40L
mkXY <- function(seed) {
  set.seed(seed)
  x <- matrix(runif(n * 2L), n, 2L)
  list(
    x = x,
    y = x[, 1L] - 0.5 * x[, 2L] + rnorm(n, sd = 0.3),
    yBin = rbinom(n, 1L, plogis(x[, 1L])),
    yOrd = ordered(
      c("lo", "mid", "hi")[1L + (x[, 1L] > 1 / 3) + (x[, 1L] > 2 / 3)],
      levels = c("lo", "mid", "hi")
    ),
    yNb = rnbinom(n, size = 5L, mu = exp(x[, 1L])),
    time = exp(x[, 1L] + rnorm(n, sd = 0.3)),
    status = rep(1, n),
    z = rbinom(n, 1L, 0.5),
    group = sample.int(4L, n, replace = TRUE),
    label = sample.int(3L, n, replace = TRUE)
  )
}

ctl <- function(seed, ...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 4L,
    n.burn = 0L,
    n.samples = 1L,
    updateState = TRUE,
    seed = seed,
    ...
  )
}

varForest <- function() dbartsForests$varianceForest(n.trees = 3L)

## ---- 3. Per-variant base fixture (the R5 dbartsSampler the mutations act on)

# variants: the twelve row keys, plus the qualified sub-cases the cells name -
# aft.var (aft under a variance forest), bcf.probit and bcf.logistic (bcf over
# a latent family; plain "bcf" is the gaussian one) and hazard.logistic (the
# third hazard spelling). hurdle has no sampler: dbarts() refuses it, and that
# refusal is itself a checked claim.
buildBase <- function(variant, seed, extra = list()) {
  d <- mkXY(seed)
  test <- d$x[1:5, , drop = FALSE]
  bcfForests <- function() {
    list(
      dbartsForests$forest(),
      dbartsForests$forest(basis = ~ factor(d$z))
    )
  }
  # no 'test =' for bcf: a treatment forest refuses test predictors outright
  # (docs/design/bcf.md), and none of the cells this fixture backs need one
  args <- switch(
    variant,
    gaussian = list(d$x, d$y, test = test, control = ctl(seed)),
    student = list(
      d$x,
      d$y,
      test = test,
      family = dbarts:::student(df = 5),
      control = ctl(seed)
    ),
    probit = list(
      d$x,
      d$yBin,
      test = test,
      family = "probit",
      control = ctl(seed)
    ),
    logistic = list(
      d$x,
      d$yBin,
      test = test,
      family = "logistic",
      control = ctl(seed)
    ),
    ordinal = list(
      d$x,
      d$yOrd,
      test = test,
      family = "ordinal",
      control = ctl(seed)
    ),
    nbinom = list(
      d$x,
      d$yNb,
      test = test,
      family = dbarts:::nbinom(shape = 5),
      control = ctl(seed)
    ),
    multinom = list(
      d$x,
      factor(d$label),
      test = test,
      family = "multinomial",
      control = ctl(seed)
    ),
    aft = list(
      d$x,
      cbind(d$time, d$status),
      test = test,
      family = "aft",
      control = ctl(seed)
    ),
    aft.var = list(
      d$x,
      cbind(d$time, d$status),
      test = test,
      family = "aft",
      variance = varForest(),
      control = ctl(seed)
    ),
    hazard = list(
      d$x,
      survival::Surv(d$time, d$status),
      family = "hazard",
      control = ctl(seed)
    ),
    hazard.logistic = list(
      d$x,
      survival::Surv(d$time, d$status),
      family = "hazard.logistic",
      control = ctl(seed)
    ),
    hurdle = list(
      d$x,
      ifelse(d$yBin == 1L, exp(d$x[, 1L]) + 0.1, 0),
      family = "hurdle.lognormal",
      control = ctl(seed)
    ),
    bcf = list(d$x, d$y, forests = bcfForests(), control = ctl(seed)),
    bcf.probit = list(
      d$x,
      d$yBin,
      forests = bcfForests(),
      family = "probit",
      control = ctl(seed)
    ),
    bcf.logistic = list(
      d$x,
      d$yBin,
      forests = bcfForests(),
      family = "logistic",
      control = ctl(seed)
    ),
    hetero = list(
      d$x,
      d$y,
      test = test,
      variance = varForest(),
      control = ctl(seed)
    ),
    stop("composition-matrix: no base fixture recipe for '", variant, "'")
  )
  if (length(extra) > 0L) {
    args[names(extra)] <- extra
  }
  sampler <- do.call(dbarts, args)
  invisible(sampler$run(0L, 1L))
  sampler
}

## ---- 4. Outcome classification --------------------------------------------
## A caught condition is REFUSES when its message reads as a deliberate,
## named stop() - this codebase's own convention throughout. It is ERRORS when
## the message instead looks like a generic R/engine malfunction: a bounds,
## conformability, symbol or NA-propagation failure nothing chose to name. A
## dbarts tombstone ("unused argument 'x' passed to 'f'") is named, so it is
## not generic; R's bare "unused argument (x = ..)" is.
genericFailure <- c(
  "subscript out of bounds",
  "non-conformable",
  "missing value where",
  "NA/NaN/Inf",
  "attempt to apply non-function",
  "could not find function",
  "unused argument",
  "is missing, with no default",
  "object '.*' not found",
  "no applicable method",
  "argument .* matches multiple"
)

# the runner's own mistake (an invalid fixture or argument), never a verdict
harnessStop <- function(...) {
  stop(structure(
    class = c("cmxHarness", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# an S cell's effect did not show: reads as a disagreement, not a refusal
effect <- function(ok, what) {
  if (!isTRUE(ok)) {
    stop("effect missing: ", what, call. = FALSE)
  }
  invisible(TRUE)
}

attempt <- function(expr) {
  ok <- tryCatch(
    {
      force(expr)
      TRUE
    },
    cmxHarness = function(e) structure(conditionMessage(e), harness = TRUE),
    error = function(e) conditionMessage(e)
  )
  if (isTRUE(ok)) {
    return(list(verdict = "CONSTRUCTS", detail = ""))
  }
  if (isTRUE(attr(ok, "harness"))) {
    return(list(verdict = "HARNESS", detail = as.character(ok)))
  }
  isGeneric <- !grepl("passed to '", ok, fixed = TRUE) &&
    any(vapply(genericFailure, grepl, logical(1L), x = ok))
  list(verdict = if (isGeneric) "ERRORS" else "REFUSES", detail = ok)
}

## ---- 5. Construction doors -------------------------------------------------

doorArgs <- function(...) {
  list(
    n.trees = 3L,
    n.burn = 0L,
    n.samples = 1L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  )
}

# the bart() argument recipe shared by the bart door, extract(loglik) and the
# alternate-family arcs; NULL for bcf, whose forest() term is a formula
bartArgs <- function(family, d) {
  switch(
    family,
    gaussian = list(d$x, d$y),
    student = list(d$x, d$y, family = dbarts:::student(df = 5)),
    probit = list(d$x, d$yBin, family = "probit"),
    logistic = list(d$x, d$yBin, family = "logistic"),
    ordinal = list(d$x, d$yOrd, family = "ordinal"),
    nbinom = list(d$x, d$yNb, family = "nbinom"),
    multinom = list(d$x, factor(d$label), family = "multinomial"),
    aft = list(d$x, cbind(d$time, d$status), family = "aft"),
    hazard = list(d$x, survival::Surv(d$time, d$status), family = "hazard"),
    hurdle = list(
      d$x,
      ifelse(d$yBin == 1L, exp(d$x[, 1L]) + 0.1, 0),
      family = "hurdle.lognormal"
    ),
    hetero = list(d$x, d$y, variance = varForest()),
    NULL
  )
}

doBart <- function(family, seed, ...) {
  d <- mkXY(seed)
  if (family == "bcf") {
    df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L], z = d$z, y = d$y)
    return(do.call(
      bart,
      c(
        list(y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z)), df),
        doorArgs(...)
      )
    ))
  }
  do.call(bart, c(bartArgs(family, d), doorArgs(...)))
}

# bartBT() takes 0.9-x's argument names and has no family formal, so a family
# is passed as one anyway; a binary response is fitted as probit
doBartBT <- function(family, seed) {
  d <- mkXY(seed)
  a <- list(
    ntree = 3L,
    ndpost = 1L,
    nskip = 0L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE
  )
  y <- if (family == "probit") d$yBin else d$y
  famArg <- if (family %in% c("gaussian", "probit")) {
    list()
  } else {
    list(family = family)
  }
  do.call(bartBT, c(list(d$x, y), famArg, a))
}

doXbart <- function(family, seed, response = NULL) {
  d <- mkXY(seed)
  df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L], y = d$y, yBin = d$yBin)
  a <- list(
    n.samples = 2L,
    n.reps = 1L,
    n.burn = c(1L, 1L),
    n.trees = 3L,
    n.threads = 1L,
    verbose = FALSE
  )
  call <- function(f, ...) do.call(xbart, c(list(f, df), a, list(...)))
  if (!is.null(response)) {
    df$y <- switch(
      response,
      ordered = d$yOrd,
      factor = factor(d$label),
      harnessStop("no xbart() response recipe '", response, "'")
    )
    return(call(y ~ x1 + x2))
  }
  switch(
    family,
    gaussian = call(y ~ x1 + x2),
    probit = call(yBin ~ x1 + x2, family = "probit"),
    logistic = call(yBin ~ x1 + x2, family = "logistic"),
    student = call(y ~ x1 + x2, family = dbarts:::student(df = 5)),
    ordinal = call(y ~ x1 + x2, family = "ordinal"),
    nbinom = call(y ~ x1 + x2, family = "nbinom"),
    multinom = call(y ~ x1 + x2, family = "multinomial"),
    aft = call(y ~ x1 + x2, family = "aft"),
    hazard = call(y ~ x1 + x2, family = "hazard"),
    hurdle = call(y ~ x1 + x2, family = "hurdle.lognormal"),
    bcf = call(y ~ x1 + x2 + forest(x1 + x2, basis = ~yBin)),
    hetero = call(y ~ x1 + x2, variance = varForest()),
    harnessStop("no xbart() recipe for '", family, "'")
  )
}

# the flat C dbarts.h route: the same compiled consumer test-capi.R drives,
# resolved through R_GetCCallable exactly as a LinkingTo package would. The
# header creates no sampler, so the handle is the address in an R-built
# sampler's external pointer; the probe reads the sampler's dimensions
# (observations, predictors, fitted values and variable-count forests per
# observation) and runs it through the flat entry.
capiConsumer <- local({
  helper <- system.file("common", "capiConsumer.R", package = "dbarts")
  if (!nzchar(helper)) {
    return(NULL)
  }
  env <- new.env()
  source(helper, local = env)
  consumer <- env$compileCapiConsumer("capi", "the C API consumer")
  if (!is.null(consumer$skip)) NULL else consumer
})
if (is.null(capiConsumer)) {
  stop(
    "composition-matrix: the flat C consumer could not be built (it needs ",
    "the installed dbarts with inst/tinytest/capi and a C toolchain); the ",
    "flat C cells cannot be judged without it"
  )
}

doFlatC <- function(sampler, variant) {
  ptr <- sampler$getPointer()
  dims <- capiConsumer$CALL("capi_dims", ptr)
  effect(
    dims[1L] == NROW(sampler$data@y),
    "flat dimensions report the sampler's observations"
  )
  if (variant == "multinom") {
    # footnote f4: K fitted values per observation, one count set per forest
    effect(dims[9L] == 3L, "flat multinomial: K values per observation")
    effect(dims[10L] == 3L, "flat multinomial: K variable-count forests")
  }
  effect(
    isTRUE(capiConsumer$CALL("capi_run_plain", ptr, 0L, 1L)),
    "flat run"
  )
}

## ---- 6. Mutation channels (table 1) and universal capabilities ------------

# each takes the sampler and the variant it was built as
mutate <- list(
  setWeights = function(s, v) {
    nr <- length(s$data@y)
    w <- if (grepl("logistic", v, fixed = TRUE)) {
      sample(1:3, nr, replace = TRUE)
    } else {
      runif(nr, 0.5, 1.5)
    }
    if (length(w) != nr) {
      harnessStop("weights of length ", length(w), " for ", nr, " rows")
    }
    s$setWeights(w)
    effect(isTRUE(all.equal(s$data@weights, w)), "weights read back")
  },
  setSigma = function(s, v) {
    s$setSigma(1.7)
    effect(isTRUE(all.equal(as.numeric(s$getSigmas()), 1.7)), "sigma read back")
  },
  getLatents = function(s, v) {
    latents <- s$getLatents()
    # the matrix's multinomial cell is R by a DECIDED decline that "reports
    # nothing" (docs/design/multinomial.md): the refusal is a NULL, not a stop
    if (v == "multinom" && is.null(latents)) {
      stop("$getLatents reports NULL on a multinomial sampler")
    }
    effect(
      is.numeric(latents) && length(latents) == length(s$data@y),
      "latents are one number per row"
    )
  },
  nullLatents = function(s, v) {
    if (!is.null(s$getLatents())) {
      return(invisible(TRUE))
    }
    stop("$getLatents reports NULL, as the '-' cell says", call. = FALSE)
  },
  updateScale = function(s, v) {
    s$setOffset(rep(0, length(s$data@y)), updateScale = TRUE)
  },
  updateScaleActiveRows = function(s, v) {
    nr <- length(s$data@y)
    s$setActiveRows(c(rep(1, nr - 1L), 0), updateScale = TRUE)
  },
  updateScaleWeights = function(s, v) {
    nr <- length(s$data@y)
    w <- if (grepl("logistic", v, fixed = TRUE)) {
      sample(1:3, nr, replace = TRUE)
    } else {
      runif(nr, 0.5, 1.5)
    }
    s$setWeights(w, updateScale = TRUE)
  },
  testSurface = function(s, v) {
    x <- s$data@x
    s$setTestPredictor(x[seq_len(min(5L, nrow(x))), , drop = FALSE])
    effect(nrow(s$data@x.test) == min(5L, nrow(x)), "test predictors read back")
  },
  setResponse = function(s, v) {
    y <- s$data@y
    if (v %in% c("gaussian", "student", "hetero", "bcf")) {
      s$setResponse(y + 1)
      effect(isTRUE(all.equal(s$data@y, y + 1)), "response read back")
    } else {
      s$setResponse(y)
    }
  },
  setOffset = function(s, v) {
    s$setOffset(rep(0.25, length(s$data@y)))
    effect(all(s$data@offset == 0.25), "offset read back")
  },
  setPredictor = function(s, v) {
    x <- s$data@x
    s$setPredictor(x)
    effect(isTRUE(all.equal(s$data@x, x)), "predictors read back")
  },
  activeRows = function(s, v) {
    nr <- length(s$data@y)
    mask <- c(rep(1, nr - 1L), 0)
    s$setActiveRows(mask)
    effect(
      identical(as.numeric(s$activeRows), mask),
      "mask held by the sampler"
    )
  },
  # a named sd written over the sampler's own prior, fixed where its k is
  # fixed and under invchi() where its k is drawn; the engine states either
  # as k against the data's scale, so the probe is the named-sd writer
  leafSd = function(s, v) {
    # the couplings carry no single leaf prior to read; the writer is the probe
    if (v %in% c("bcf", "multinom")) {
      return(s$setLeafPrior(dbarts:::normal(sd = 0.5)))
    }
    before <- s$getLeafPrior()$leaf.prior@k
    if (is.numeric(before)) {
      s$setLeafPrior(dbarts:::normal(sd = 1.2 / before))
    } else {
      s$setLeafPrior(dbarts:::normal(
        sd = dbarts:::invchi(before@degreesOfFreedom, 0.6)
      ))
    }
    effect(
      !identical(s$getLeafPrior()$leaf.prior@k, before),
      "leaf prior changed by the named sd"
    )
  },
  growFromRoot = function(s, v) s$growFromRoot(1L),
  setData = function(s, v) {
    if (v == "gaussian") {
      s$setData(dbartsData(s$data@x[1:30, ], s$data@y[1:30]))
      effect(length(s$data@y) == 30L, "whole-data replacement read back")
    } else {
      s$setData(s$data)
    }
  },
  setForestWeights = function(s, v) {
    s$setForestWeights(1L, c(rep(1, length(s$data@y) - 1L), 0))
  }
)

# how a table 1 or 2 qualifier picks fixtures; "" is the row itself (the
# gaussian bcf for plain "bcf"). A qualifier absent here stops the run.
qualifierVariants <- list(
  "bcf|setWeights|gaussian/logistic" = c("bcf", "bcf.logistic"),
  "bcf|setSigma|gaussian only" = "bcf",
  "bcf|getLatents|probit/logistic only" = c("bcf.probit", "bcf.logistic"),
  "aft|setSigma|no variance forest" = "aft",
  "aft|setSigma|variance forest" = "aft.var",
  "aft|updateScale|a variance forest re-anchored with it" = c("aft", "aft.var"),
  "hetero|setWeights|gaussian" = "hetero",
  "hetero|setWeights|aft" = "aft.var",
  "hetero|getLatents|gaussian" = "hetero",
  "hetero|getLatents|aft" = "aft.var",
  "hazard|flatc|as probit" = c("hazard", "hazard.logistic")
)

# the package's own refusal text for every R or M cell, matched as a regular
# expression against the message; keyed "row|capability" or, where the cell is
# qualified, "row|capability|qualifier". A cell with no entry stops the run.
multinomText <- "is not available on a multinomial sampler"
caseWeightsText <- "models do not support case weights other than 0 and 1"
fixedSigmaText <- "fixes the residual standard deviation by definition"
varOwnsText <- "variance forest owns the residual scale"
tableText <- c(
  "probit|setWeights" = caseWeightsText,
  "ordinal|setWeights" = caseWeightsText,
  "nbinom|setWeights" = caseWeightsText,
  "hazard|setWeights" = caseWeightsText,
  "aft|setWeights" = "aft \\(survival\\) models do not support case weights",
  "hetero|setWeights|aft" = "aft \\(survival\\) models do not support case weights",
  "multinom|setWeights" = multinomText,
  "probit|setSigma" = fixedSigmaText,
  "logistic|setSigma" = fixedSigmaText,
  "ordinal|setSigma" = fixedSigmaText,
  "nbinom|setSigma" = fixedSigmaText,
  "hazard|setSigma" = fixedSigmaText,
  "aft|setSigma|variance forest" = varOwnsText,
  "hetero|setSigma" = varOwnsText,
  "multinom|setSigma" = multinomText,
  "multinom|getLatents" = "reports NULL on a multinomial sampler",
  "multinom|updateScale" = multinomText,
  "bcf|updateScale" = "does not support a sampler that carries forest amplitudes",
  "bcf|testSurface" = "combined test fits are undefined",
  "student|xbart" = "does not fit family",
  "ordinal|xbart" = "'family' should be one of",
  "nbinom|xbart" = "'family' should be one of",
  "multinom|xbart" = "'family' should be one of",
  "aft|xbart" = "'family' should be one of",
  "hazard|xbart" = "'family' should be one of",
  "hurdle|xbart" = "'family' should be one of",
  "bcf|xbart" = "does not take a forest\\(\\) term",
  "hetero|xbart" = "unused argument 'variance' passed to 'xbart'",
  "hurdle|flatc" = "composes two independent samplers"
)
lookupText <- function(fam, cap, qualifier) {
  keys <- c(paste(fam, cap, qualifier, sep = "|"), paste(fam, cap, sep = "|"))
  hit <- keys[keys %in% names(tableText)]
  if (length(hit) == 0L) {
    stop(
      "composition-matrix: no refusal text pinned for the cell '",
      keys[1L],
      "' - matrix vocabulary needs a human look"
    )
  }
  tableText[[hit[1L]]]
}

## ---- 7. The checks ---------------------------------------------------------

checks <- list()
# act(setupResult, seed) runs inside attempt(); setup(seed) does not, so a
# fixture that fails to build is a harness failure rather than a refusal -
# unless the claim is S (then the failure is the disagreement) or the check
# says the fixture's construction is itself the claim (setupClaim). msg is the
# regular expression an R or M check's refusal must match; generic lets a
# bare R error carry that text (a cell unreachable for want of a formal).
addCheck <- function(
  source,
  label,
  expect,
  act,
  setup = NULL,
  anchor = NULL,
  setupClaim = FALSE,
  msg = NULL,
  generic = FALSE
) {
  if (expect != "S" && is.null(msg)) {
    stop("composition-matrix: check '", label, "' needs a refusal text")
  }
  checks[[length(checks) + 1L]] <<- list(
    source = source,
    label = label,
    expect = expect,
    act = act,
    setup = setup,
    anchor = anchor,
    setupClaim = setupClaim,
    msg = msg,
    generic = generic
  )
}
addMutation <- function(
  source,
  label,
  expect,
  variant,
  probe,
  anchor = NULL,
  msg = NULL
) {
  force(variant)
  force(probe)
  addCheck(
    source,
    label,
    expect,
    act = function(s, seed) mutate[[probe]](s, variant),
    setup = function(seed) buildBase(variant, seed),
    anchor = anchor,
    msg = msg
  )
}

lines <- readLines(matrixPath)
matrixText <- gsub("\\s+", " ", paste(lines, collapse = " "))

# -- tables 1 and 2, read from the file
invisible(lapply(parseTables(lines), function(cl) {
  fam <- cl$family
  cap <- cl$capability
  label <- sprintf("%-8s %-11s %s", fam, cap, cl$token)
  if (cl$qualifier != "") {
    label <- paste0(label, ", ", cl$qualifier)
  }
  if (fam == "hurdle" && cl$token == "-") {
    return(invisible()) # no sampler, no concept: the row's footnote [f10]
  }
  if (cl$token == "-") {
    if (cap == "updateScale") {
      # [f7]: ignored rather than refused where the transform is fixed
      addMutation(cl$table, label, "S", fam, cap)
    } else if (cap == "getLatents" && fam %in% c("gaussian", "hetero")) {
      # reports nothing: the latent concept does not apply to this sampler
      addMutation(
        cl$table,
        label,
        "R",
        fam,
        "nullLatents",
        msg = "reports NULL, as the '-' cell says"
      )
    }
    return(invisible())
  }
  variants <- if (cl$qualifier == "") {
    fam
  } else {
    key <- paste(fam, cap, cl$qualifier, sep = "|")
    if (is.null(qualifierVariants[[key]])) {
      stop(
        "composition-matrix: no fixture for qualifier '",
        key,
        "' - matrix vocabulary needs a human look"
      )
    }
    qualifierVariants[[key]]
  }
  text <- if (cl$token == "S") NULL else lookupText(fam, cap, cl$qualifier)
  for (v in variants) {
    local({
      v <- v
      lab <- if (length(variants) > 1L) paste0(label, " [", v, "]") else label
      if (cap == "xbart") {
        addCheck(
          cl$table,
          lab,
          cl$token,
          act = function(s, seed) doXbart(fam, seed),
          msg = text
        )
      } else if (cap == "flatc") {
        addCheck(
          cl$table,
          lab,
          cl$token,
          act = function(s, seed) doFlatC(s, v),
          setup = function(seed) buildBase(v, seed),
          setupClaim = fam == "hurdle",
          msg = text
        )
      } else {
        addMutation(cl$table, lab, cl$token, v, cap, msg = text)
      }
    })
  }
}))

# xbart's other refusal of the classification families: an ordered or
# categorical response, whose message redirects to bart()
for (r in c("ordered", "factor")) {
  local({
    r <- r
    addCheck(
      "table2",
      sprintf(
        "%-8s xbart       R, %s response",
        if (r == "ordered") "ordinal" else "multinom",
        r
      ),
      "R",
      act = function(s, seed) doXbart("gaussian", seed, response = r),
      msg = "xbart does not fit a 3-level"
    )
  })
}

# -- prose: construction reach
reachAnchor <- "gaussian and probit reach all three"
onlyAnchor <- "student, logistic and aft reach `bart()` and `dbarts()` only"
bartBTAnchor <- "are out of `bartBT()`'s reach"
bartBTText <- "unused argument \\(family = "
for (fam in c("gaussian", "probit")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s bartBT", fam),
      "S",
      function(s, seed) doBartBT(fam, seed),
      anchor = reachAnchor
    )
  })
}
for (fam in c("gaussian", "student", "probit", "logistic", "aft")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s bart", fam),
      "S",
      function(s, seed) doBart(fam, seed),
      anchor = reachAnchor
    )
    addCheck(
      "prose",
      sprintf("reach %-8s dbarts", fam),
      "S",
      function(s, seed) NULL,
      setup = function(seed) buildBase(fam, seed),
      anchor = reachAnchor
    )
  })
}
for (fam in c("student", "logistic", "aft")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s bartBT", fam),
      "M",
      function(s, seed) doBartBT(fam, seed),
      anchor = onlyAnchor,
      msg = bartBTText,
      generic = TRUE
    )
  })
}
for (fam in c("ordinal", "nbinom", "multinom", "hurdle", "hazard", "bcf")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s bartBT", fam),
      "M",
      function(s, seed) doBartBT(fam, seed),
      anchor = bartBTAnchor,
      msg = bartBTText,
      generic = TRUE
    )
  })
}
addCheck(
  "prose",
  "reach hetero   bartBT",
  "M",
  function(s, seed) {
    d <- mkXY(seed)
    bartBT(d$x, d$y, variance = varForest())
  },
  anchor = bartBTAnchor,
  msg = "unused argument \\(variance = ",
  generic = TRUE
)
for (fam in c(
  "ordinal",
  "nbinom",
  "multinom",
  "hurdle",
  "hazard",
  "bcf",
  "hetero"
)) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s bart", fam),
      "S",
      function(s, seed) doBart(fam, seed),
      anchor = "`bart()` ships all six"
    )
  })
}
for (fam in c("ordinal", "nbinom", "multinom", "hazard", "bcf", "hetero")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("reach %-8s dbarts", fam),
      "S",
      function(s, seed) NULL,
      setup = function(seed) buildBase(fam, seed),
      anchor = "`dbarts()` + R5 ships every one but hurdle"
    )
  })
}
addCheck(
  "prose",
  "reach hurdle   dbarts",
  "R",
  function(s, seed) buildBase("hurdle", seed),
  anchor = "`dbarts()` + R5 ships every one but hurdle",
  msg = "composes two independent samplers"
)
addCheck(
  "prose",
  "reach twopart  bart",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    do.call(bart, c(list(d$x, d$y, family = "twopart"), doorArgs()))
  },
  anchor = "is not a family token at either door",
  msg = "'family' should be one of"
)
addCheck(
  "prose",
  "bartBT factor of three levels",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    bartBT(d$x, factor(d$label), ntree = 3L, ndpost = 1L, nskip = 0L)
  },
  anchor = "a factor of three or more levels",
  msg = "does not fit a factor response with three or more levels"
)

# -- prose: forests = admits gaussian/probit/logistic only; variance forest
# requires gaussian or aft
for (fam in c("aft", "ordinal", "nbinom")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("forests= %-8s dbarts", fam),
      "R",
      function(s, seed) {
        d <- mkXY(seed)
        y <- switch(
          fam,
          aft = cbind(d$time, d$status),
          ordinal = d$yOrd,
          d$yNb
        )
        dbarts(
          d$x,
          y,
          forests = list(
            dbartsForests$forest(),
            dbartsForests$forest(basis = ~ factor(d$z))
          ),
          family = fam,
          control = ctl(seed)
        )
      },
      anchor = "A `forests =` fit resolves gaussian, probit or logistic only",
      msg = paste0("a treatment forest does not support family \"", fam, "\"")
    )
  })
}
addCheck(
  "prose",
  "forests= aft      bart forest() term",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L], z = d$z, y = d$time)
    do.call(
      bart,
      c(
        list(
          survival::Surv(y, rep(1, n)) ~ x1 + forest(x1, basis = ~z),
          df,
          family = "aft"
        ),
        doorArgs()
      )
    )
  },
  anchor = "under an identical gate",
  msg = "does not support a forest\\(\\) formula term"
)
# the string "extra:variance" below is cited by docs/design/aft-variance-forest.md
for (fam in c("gaussian", "aft")) {
  local({
    fam <- fam
    v <- if (fam == "aft") "aft.var" else "hetero"
    addCheck(
      "prose",
      sprintf("variance= %-8s", fam),
      "S",
      function(s, seed) NULL,
      setup = function(seed) buildBase(v, seed),
      anchor = "requires `family = \"gaussian\"` or `\"aft\"`"
    )
  })
}
for (fam in c("probit", "logistic", "ordinal", "nbinom", "student")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("variance= %-8s", fam),
      "R",
      function(s, seed) buildBase(fam, seed, list(variance = varForest())),
      anchor = if (fam == "student") {
        "it refuses Student-t residuals"
      } else {
        "The other four single-forest families"
      },
      msg = if (fam == "student") {
        "a variance forest does not support Student-t residuals"
      } else {
        paste0(
          "a variance forest requires family = \"gaussian\" or \"aft\"; ",
          "family \"",
          fam,
          "\""
        )
      }
    )
  })
}
addCheck(
  "prose",
  "variance= gaussian monotone",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    do.call(
      bart,
      c(
        list(d$x, d$y, variance = varForest(), monotone = c(1L, 0L)),
        doorArgs()
      )
    )
  },
  anchor = "not supported with monotone constraints",
  msg = "a variance forest is not supported with monotone constraints"
)
# the variance forest never takes DART itself; a DART mean forest composes
addCheck(
  "prose",
  "variance= gaussian dart mean forest",
  "S",
  function(s, seed) NULL,
  setup = function(seed) {
    buildBase("hetero", seed, list(split.probs = dbarts:::dart()))
  },
  anchor = "It never takes DART either"
)

# -- prose: DART refused for both couplings
dartAnchor <- "DART is refused for both multi-forest couplings by value"
addCheck(
  "prose",
  "dart bcf",
  "R",
  function(s, seed) {
    buildBase("bcf", seed, list(split.probs = dbarts:::dart()))
  },
  anchor = dartAnchor,
  msg = "a model of several forests does not support 'split.probs'"
)
addCheck(
  "prose",
  "dart multinom",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    do.call(
      bart,
      c(bartArgs("multinom", d), doorArgs(split.probs = dbarts:::dart()))
    )
  },
  anchor = dartAnchor,
  msg = "multinomial .softmax. model does not support .split.probs."
)

# -- prose: warm start and grow-from-root
warmAnchor <- paste0(
  "Warm start and grow-from-root are refused, as `M` rather than `R`, ",
  "for the four alternate-family `bart` arcs"
)
for (fam in c("ordinal", "nbinom", "multinom", "hurdle")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("grow-from-root %-8s bart", fam),
      "M",
      function(s, seed) {
        d <- mkXY(seed)
        do.call(bart, c(bartArgs(fam, d), doorArgs(n.grow.sweeps = 1L)))
      },
      anchor = warmAnchor,
      msg = "does not support 'warm.start' or 'n.grow.sweeps'"
    )
  })
}
treeKey <- c("tree", "var", "value")
for (fam in c("gaussian", "student", "probit", "logistic", "aft", "hazard")) {
  local({
    fam <- fam
    addMutation(
      "prose",
      sprintf("grow-from-root %-8s", fam),
      "S",
      fam,
      "growFromRoot",
      anchor = "grow-from-root itself ships"
    )
    addCheck(
      "prose",
      sprintf("warm start %-8s", fam),
      "S",
      function(s, seed) {
        s$sampler$installTrees(s$donor)
        effect(
          identical(
            s$sampler$getTrees(current = TRUE)[treeKey],
            s$donor$getTrees(current = TRUE)[treeKey]
          ),
          "installed trees equal the donor's"
        )
      },
      setup = function(seed) {
        list(
          sampler = buildBase(fam, seed),
          # same data, another stream: the donor's trees differ from the sampler's
          donor = buildBase(fam, seed, list(control = ctl(seed + 1000L)))
        )
      },
      anchor = warmAnchor
    )
  })
}
addMutation(
  "prose",
  "grow-from-root bcf (two forests)",
  "S",
  "bcf",
  "growFromRoot",
  anchor = "covered at two forests"
)
addMutation(
  "prose",
  "grow-from-root hetero",
  "S",
  "hetero",
  "growFromRoot",
  anchor = "grow-from-root itself ships"
)
addCheck(
  "prose",
  "warm start bcf donor (multi-forest)",
  "M",
  function(s, seed) s$sampler$installTrees(s$donor),
  setup = function(seed) {
    list(
      sampler = buildBase("bcf", seed),
      donor = buildBase("bcf", seed + 1000L)
    )
  },
  anchor = "and for the multi-forest donor warm start",
  msg = "does not support a multi-forest sampler"
)

# -- prose: setResponse / setOffset / setPredictor open wherever a sampler
# exists, multinomial blocking the first two
openAnchor <- "are open wherever a sampler exists"
samplerVariants <- c(
  "gaussian",
  "student",
  "probit",
  "logistic",
  "ordinal",
  "nbinom",
  "multinom",
  "aft",
  "hazard",
  "bcf",
  "hetero"
)
for (v in samplerVariants) {
  for (probe in c("setResponse", "setOffset", "setPredictor")) {
    blocked <- v == "multinom" && probe != "setPredictor"
    addMutation(
      "prose",
      sprintf("%-8s %s", v, probe),
      if (blocked) "R" else "S",
      v,
      probe,
      anchor = openAnchor,
      msg = if (blocked) multinomText else NULL
    )
  }
}

# -- prose: updateScale = TRUE on setActiveRows (every variant) and setWeights
# (where the case-weights cell is S) follows the unit-scale cell: run where it
# is a no-op or a re-anchor, refused on multinomial and bcf
scaleAnchor <- paste0(
  "`updateScale = TRUE`, on `setResponse`, `setOffset`, `setWeights` and ",
  "`setActiveRows`, follows the unit-scale column"
)
for (v in samplerVariants) {
  refused <- v %in% c("multinom", "bcf")
  for (probe in c("updateScaleActiveRows", "updateScaleWeights")) {
    if (
      probe == "updateScaleWeights" &&
        !v %in% c("gaussian", "student", "logistic", "hetero", "bcf")
    ) {
      next
    }
    addMutation(
      "prose",
      sprintf("%-8s %s", v, probe),
      if (refused) "R" else "S",
      v,
      probe,
      anchor = scaleAnchor,
      msg = if (v == "multinom") {
        multinomText
      } else if (refused) {
        "does not support a sampler that carries forest amplitudes"
      } else {
        NULL
      }
    )
  }
}

# -- prose: whole-data setData needs a sampler that owns its predictors, one
# forest, and a status-free family
dataAnchor <- "Whole-data `setData` (n free) is narrower than `setPredictor`"
addMutation(
  "prose",
  "setData gaussian",
  "S",
  "gaussian",
  "setData",
  anchor = dataAnchor
)
addMutation(
  "prose",
  "setData aft",
  "R",
  "aft",
  "setData",
  anchor = dataAnchor,
  msg = "fix the censoring structure at creation"
)
addMutation(
  "prose",
  "setData bcf",
  "R",
  "bcf",
  "setData",
  anchor = dataAnchor,
  msg = "does not support a sampler that carries forest amplitudes"
)
addMutation(
  "prose",
  "setData multinom",
  "R",
  "multinom",
  "setData",
  anchor = dataAnchor,
  msg = multinomText
)

# -- prose: bcf "only" cells, the third hazard spelling
addMutation(
  "prose",
  "bcf.probit setSigma",
  "R",
  "bcf.probit",
  "setSigma",
  anchor = "S, gaussian only",
  msg = fixedSigmaText
)
addMutation(
  "prose",
  "bcf.logistic setSigma",
  "R",
  "bcf.logistic",
  "setSigma",
  anchor = "S, gaussian only",
  msg = fixedSigmaText
)
addMutation(
  "prose",
  "bcf getLatents gaussian",
  "R",
  "bcf",
  "nullLatents",
  anchor = "S, probit/logistic only",
  msg = "reports NULL, as the '-' cell says"
)
addMutation(
  "prose",
  "hazard.logistic setWeights",
  "S",
  "hazard.logistic",
  "setWeights",
  anchor = "or, under the third spelling, the logistic row"
)
addMutation(
  "prose",
  "hazard.logistic getLatents",
  "S",
  "hazard.logistic",
  "getLatents",
  anchor = "or, under the third spelling, the logistic row"
)

# -- prose: zero weights follow the case-weights column except logistic;
# probit, ordinal and nbinom take a 0/1 weight vector at creation
zeroAnchor <- "A zero weight, which takes a row out of the likelihood"
for (v in c(
  "gaussian",
  "student",
  "probit",
  "ordinal",
  "hetero",
  "bcf",
  "logistic",
  "bcf.logistic",
  "aft",
  "nbinom"
)) {
  local({
    v <- v
    refused <- v %in% c("logistic", "bcf.logistic", "aft")
    addCheck(
      "prose",
      sprintf("zero weight at creation %-8s", v),
      if (refused) "R" else "S",
      function(s, seed) {
        buildBase(v, seed, list(weights = c(0, rep(1, n - 1L))))
      },
      anchor = zeroAnchor,
      msg = if (grepl("logistic", v, fixed = TRUE)) {
        "logistic weights are observation counts and must be positive integers"
      } else {
        "models do not support case weights"
      }
    )
  })
}

# -- prose: the universal capabilities
for (v in c(samplerVariants, "aft.var", "bcf.probit", "bcf.logistic")) {
  addMutation(
    "prose",
    sprintf("active-rows mask %-12s", v),
    "S",
    v,
    "activeRows",
    anchor = "the mid-chain active-rows mask"
  )
}
addMutation(
  "prose",
  "per-forest mask multinom",
  "R",
  "multinom",
  "setForestWeights",
  anchor = "Multinomial's per-forest active-rows mask: refused permanently",
  msg = multinomText
)
for (fam in c(
  "gaussian",
  "student",
  "probit",
  "logistic",
  "ordinal",
  "nbinom",
  "multinom",
  "aft",
  "hazard",
  "hurdle",
  "bcf",
  "hetero"
)) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("extract loglik %-8s", fam),
      "S",
      function(s, seed) {
        ll <- extract(s, type = "loglik")
        effect(is.numeric(ll) && length(ll) > 0L, "loglik is numeric")
      },
      setup = function(seed) doBart(fam, seed),
      anchor = "`extract(type = \"loglik\")` (every row"
    )
  })
}
for (fam in c("weights", "subset", "offset")) {
  local({
    fam <- fam
    addCheck(
      "prose",
      sprintf("hurdle bart %s", fam),
      "R",
      function(s, seed) {
        d <- mkXY(seed)
        extra <- switch(
          fam,
          weights = list(weights = rep(1, n)),
          subset = list(subset = 1:30),
          offset = list(offset = rep(0, n))
        )
        do.call(bart, c(bartArgs("hurdle", d), doorArgs(), extra))
      },
      anchor = "does not support 'weights'",
      msg = if (fam == "offset") {
        "does not support 'offset'/'offset.test'"
      } else {
        paste0("does not support '", fam, "'")
      }
    )
  })
}
for (v in c(
  "gaussian",
  "student",
  "probit",
  "logistic",
  "ordinal",
  "nbinom",
  "aft",
  "hazard",
  "hetero",
  "bcf",
  "multinom"
)) {
  coupled <- v %in% c("bcf", "multinom")
  addMutation(
    "prose",
    sprintf("named leaf sd %-8s", v),
    if (coupled) "R" else "S",
    v,
    "leafSd",
    anchor = "A named leaf-prior sd",
    msg = if (v == "bcf") {
      "takes forests = : the multi-forest calibration map"
    } else if (coupled) {
      "takes normal\\(k = \\) with a fixed k"
    } else {
      NULL
    }
  )
}

for (ck in checks) {
  if (!is.null(ck$anchor) && !grepl(ck$anchor, matrixText, fixed = TRUE)) {
    stop(
      "composition-matrix: the matrix no longer says '",
      ck$anchor,
      "' (check '",
      ck$label,
      "') - the prose probe needs a human look"
    )
  }
}
labels <- vapply(checks, function(ck) ck$label, "")
if (anyDuplicated(labels)) {
  stop(
    "composition-matrix: duplicate check label '",
    labels[anyDuplicated(labels)],
    "'"
  )
}
perPart <- table(factor(
  vapply(checks, function(ck) ck$source, ""),
  levels = names(pinChecks)
))
for (part in names(pinChecks)) {
  if (perPart[[part]] != pinChecks[[part]]) {
    pinStop(
      paste0("the number of ", part, " checks"),
      perPart[[part]],
      pinChecks[[part]]
    )
  }
}

## ---- 8. Run every check, diff against the matrix --------------------------

cat(sprintf(
  "composition-matrix: %d checks (%d table 1, %d table 2, %d prose)\n",
  length(checks),
  perPart[["table1"]],
  perPart[["table2"]],
  perPart[["prose"]]
))

confirmations <- 0L
notes <- character()
disagreements <- character()
harnessFailures <- character()

for (ck in checks) {
  seed <- sum(utf8ToInt(ck$label)) %% 100000L + 1L
  setupResult <- NULL
  if (!is.null(ck$setup)) {
    setupResult <- tryCatch(ck$setup(seed), error = function(e) e)
    if (inherits(setupResult, "error")) {
      if (ck$expect != "S" && !ck$setupClaim) {
        harnessFailures <- c(
          harnessFailures,
          sprintf(
            "[%s] %s: fixture did not build (%s)",
            ck$source,
            ck$label,
            conditionMessage(setupResult)
          )
        )
        next
      }
      result <- attempt(stop(conditionMessage(setupResult)))
    } else {
      result <- attempt(ck$act(setupResult, seed))
    }
  } else {
    result <- attempt(ck$act(NULL, seed))
  }
  if (result$verdict == "HARNESS") {
    harnessFailures <- c(
      harnessFailures,
      sprintf("[%s] %s: %s", ck$source, ck$label, result$detail)
    )
    next
  }
  expect <- ck$expect
  msg <- ck$msg
  if (nzchar(flipText) && grepl(flipText, ck$label, fixed = TRUE)) {
    expect <- if (expect == "S") "R" else "S"
    msg <- ""
  }
  textOk <- !is.null(msg) && grepl(msg, result$detail)
  agree <- switch(
    expect,
    S = result$verdict == "CONSTRUCTS",
    R = textOk &&
      (result$verdict == "REFUSES" ||
        (ck$generic && result$verdict == "ERRORS")),
    M = textOk &&
      (result$verdict == "REFUSES" ||
        (ck$generic && result$verdict == "ERRORS"))
  )
  line <- sprintf(
    "[%s] %s: matrix %s, got %s%s",
    ck$source,
    ck$label,
    expect,
    result$verdict,
    if (nzchar(result$detail)) sprintf(" (%s)", result$detail) else ""
  )
  if (!agree && expect != "S" && result$verdict %in% c("REFUSES", "ERRORS")) {
    line <- paste0(line, sprintf(" - expected text /%s/", msg))
  }
  if (verbose) {
    cat(if (agree) "ok   " else "FAIL ", line, "\n")
  }
  if (agree) {
    confirmations <- confirmations + 1L
    if (result$verdict == "ERRORS") {
      notes <- c(notes, line)
    }
  } else {
    disagreements <- c(disagreements, line)
  }
}

cat(sprintf("\nconfirmations: %d of %d\n", confirmations, length(checks)))
cat(sprintf(
  "\nnotes, bare R error on a cell unreachable for want of a formal (%d):\n",
  length(notes)
))
for (line in notes) {
  cat("  ", line, "\n")
}
cat(sprintf("\nDISAGREEMENTS (%d):\n", length(disagreements)))
for (line in disagreements) {
  cat("  ", line, "\n")
}
cat(sprintf("\nHARNESS FAILURES (%d):\n", length(harnessFailures)))
for (line in harnessFailures) {
  cat("  ", line, "\n")
}

quit(
  status = if (length(disagreements) + length(harnessFailures) > 0L) 1L else 0L
)
