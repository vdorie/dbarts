#!/usr/bin/env Rscript

# Executable composition matrix: an oracle for docs/design/feature-matrix.md.
# Every claim the matrix makes is turned into one small composition that is
# actually attempted, and the outcome is classified as CONSTRUCTS (runs without
# error), REFUSES (a named validation stop()) or ERRORS (anything else - the
# finding class: a generic R/engine failure rather than a deliberate refusal).
# The outcome is then diffed against what the matrix says:
#
#   S  (shipped)   must CONSTRUCT
#   R  (refused)   must REFUSE by name; a CONSTRUCTS or an ERRORS disagrees
#   M  (missing)   must not CONSTRUCT; a named refusal agrees, a generic
#                  error is counted as agreeing but listed as a note
#   -  (n/a)       unchecked, except the unit-scale column, whose footnote [f7]
#                  says updateScale = TRUE is ignored rather than refused there
#
# Two kinds of claim are read. The two tables (section 1, five capability
# columns; section 2, the xbart() and flat C columns) are parsed out of the
# file cell by cell, qualifiers included ("S, gaussian/logistic", "S, no
# variance forest / R, variance forest", "as probit"), and a cell or qualifier
# this runner has no fixture for stops the run. The prose rules (construction
# reach, the universal capabilities, zero weights, variance forest, DART, warm
# start and grow-from-root) have no table, so each is a hand-written probe
# anchored to a phrase that must still be in the matrix; a phrase that has
# moved or been reworded stops the run, so the matrix cannot drift past the
# runner unnoticed.
#
# Lives in benchmarks/R/ rather than tools/: it is a local correctness gate
# against a design document, the same role equivalence.R and bench-sampler.R
# play against a baseline, not a build-freshness check. It runs in no
# workflow; run it by hand at the release candidate.
#
# Every fixture is tiny (n = 40, a handful of trees, 0 burn-in, 1 sample,
# n.threads = 1L throughout) and reseeded per check, so a run is deterministic
# and fast; nothing here is a statistical check.
#
# Usage: Rscript benchmarks/R/composition-matrix.R [--verbose] [--flip=<text>]
# Exit status is nonzero iff at least one disagreement or harness failure was
# recorded. --verbose prints every check. --flip=<text> inverts the expectation of every check whose label
# contains <text>, to show that a check can fail.

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
# transcribed once from the doc's table headers; the first element of each
# entry is the header line's start, used to find the table, the rest are the
# dispatch keys of its columns - the cell VALUES are read out of the file
tables <- list(
  structure(
    c(
      setWeights = "case weights",
      setSigma = "sigma",
      getLatents = "latents",
      updateScale = "unit-scale transform",
      testSurface = "test fits defined"
    ),
    start = "| model | case weights |"
  ),
  structure(
    c(xbart = "`xbart()`", flatc = "flat C `dbarts.h`"),
    start = "| model | `xbart()` |"
  )
)

# a cell is one or more alternatives split by " / "; each is a token (S, R, M,
# -) with an optional qualifier after a comma, or "as <row>", which copies
# that row's cell
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
  for (cols in tables) {
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
    body <- lines[seq(hdr + 2L, length.out = length(rowKeys))]
    parsed <- list()
    for (i in seq_along(rowKeys)) {
      raw <- trimws(strsplit(body[i], "|", fixed = TRUE)[[1L]][-1L])
      if (length(raw) != length(cols) + 1L || raw[1L] != rowKeys[i]) {
        stop(
          "composition-matrix: row ",
          i,
          " under '",
          start,
          "' reads '",
          raw[1L],
          "' with ",
          length(raw),
          " cells, expected '",
          rowKeys[i],
          "' with ",
          length(cols) + 1L,
          " - row drift"
        )
      }
      parsed[[rowKeys[i]]] <- lapply(raw[-1L], parseCell)
    }
    for (family in rowKeys) {
      for (j in seq_along(cols)) {
        alts <- parsed[[family]][[j]]
        if (alts[[1L]]$token == "as") {
          alts <- parsed[[alts[[1L]]$qualifier]][[j]]
        }
        for (alt in alts) {
          claims[[length(claims) + 1L]] <- list(
            family = family,
            capability = names(cols)[j],
            token = alt$token,
            qualifier = alt$qualifier
          )
        }
      }
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
# a latent family; plain "bcf" is the gaussian one). hurdle has no sampler:
# dbarts() refuses it, and that refusal is itself a checked claim.
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
## named stop() - this codebase's own convention throughout (see the matrix's
## anchors). It is ERRORS - the finding class - when the message instead looks
## like a generic R/engine malfunction: a bounds, conformability, symbol or
## NA-propagation failure nothing chose to name.
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

attempt <- function(expr) {
  ok <- tryCatch(
    {
      force(expr)
      TRUE
    },
    error = function(e) conditionMessage(e)
  )
  if (isTRUE(ok)) {
    return(list(verdict = "CONSTRUCTS", detail = ""))
  }
  isGeneric <- any(vapply(genericFailure, grepl, logical(1L), x = ok))
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

# the bart() argument recipe shared by the bart door and extract(loglik) for
# every family that reaches bart() plainly; NULL for bcf, whose forest() term
# is a formula
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

bcfBartFit <- function(d, ...) {
  df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L], z = d$z, y = d$y)
  do.call(
    bart,
    c(list(y ~ x1 + x2 + forest(x1 + x2, basis = ~z), df), doorArgs(...))
  )
}

doBart <- function(family, seed, ...) {
  d <- mkXY(seed)
  if (family == "bcf") {
    return(bcfBartFit(d, ...))
  }
  do.call(bart, c(bartArgs(family, d), doorArgs(...)))
}

# bartBT() takes 0.9-x's argument names and has no family formal, so a
# family is passed as one anyway - the matrix's reach sentence says some of
# these reach it
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
  y <- if (family %in% c("probit", "logistic")) d$yBin else d$y
  famArg <- if (family %in% c("gaussian", "probit")) {
    list()
  } else {
    list(family = family)
  }
  do.call(bartBT, c(list(d$x, y), famArg, a))
}

doXbart <- function(family, seed) {
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
    stop("composition-matrix: no xbart() recipe for '", family, "'")
  )
}

# the flat C dbarts.h route: the same compiled consumer test-capi.R drives,
# resolved through R_GetCCallable exactly as a LinkingTo package would. The
# header creates no sampler, so the handle is the address in an R-built
# sampler's external pointer; the probe reads the sampler's dimensions and
# runs it through the flat entry.
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

doFlatC <- function(sampler) {
  if (is.null(capiConsumer)) {
    stop("could not compile the flat C API consumer")
  }
  ptr <- sampler$getPointer()
  capiConsumer$CALL("capi_dims", ptr)
  capiConsumer$CALL("capi_run_plain", ptr, 0L, 1L)
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
    s$setWeights(w)
  },
  setSigma = function(s, v) s$setSigma(1.0),
  getLatents = function(s, v) {
    latents <- s$getLatents()
    # the matrix's multinomial cell is R by a DECIDED decline that "reports
    # nothing" (docs/design/multinomial.md): the refusal is a NULL, not a stop
    if (v == "multinom" && is.null(latents)) {
      stop("$getLatents reports NULL on a multinomial sampler")
    }
    latents
  },
  updateScale = function(s, v) {
    s$setOffset(rep(0, length(s$data@y)), updateScale = TRUE)
  },
  testSurface = function(s, v) {
    x <- s$data@x
    s$setTestPredictor(x[seq_len(min(5L, nrow(x))), , drop = FALSE])
  },
  setResponse = function(s, v) s$setResponse(s$data@y),
  setOffset = function(s, v) s$setOffset(rep(0, length(s$data@y))),
  setPredictor = function(s, v) s$setPredictor(s$data@x),
  activeRows = function(s, v) {
    nr <- length(s$data@y)
    s$setActiveRows(c(rep(1, nr - 1L), 0))
  },
  # a named sd written over the sampler's own prior, fixed where its k is
  # fixed and under invchi() where its k is drawn; the engine states either
  # as k against the data's scale, so the probe is the named-sd writer
  leafSd = function(s, v) {
    # the couplings carry no single leaf prior to read; the writer is the probe
    if (v %in% c("bcf", "multinom")) {
      return(s$setLeafPrior(dbarts:::normal(sd = 0.5)))
    }
    kPrior <- s$getLeafPrior()$leaf.prior@k
    if (is.numeric(kPrior)) {
      s$setLeafPrior(dbarts:::normal(sd = 1.2 / kPrior))
    } else {
      s$setLeafPrior(dbarts:::normal(
        sd = dbarts:::invchi(kPrior@degreesOfFreedom, 0.6)
      ))
    }
  },
  growFromRoot = function(s, v) s$growFromRoot(1L)
)

# how a table 1 qualifier picks fixtures; "" is the row itself (the gaussian
# bcf for plain "bcf"). A qualifier absent here stops the run.
qualifierVariants <- list(
  "bcf|setWeights|gaussian/logistic" = c("bcf", "bcf.logistic"),
  "bcf|setSigma|gaussian only" = "bcf",
  "bcf|getLatents|probit/logistic only" = c("bcf.probit", "bcf.logistic"),
  "aft|setSigma|no variance forest" = "aft",
  "aft|setSigma|variance forest" = "aft.var",
  "aft|updateScale|a variance forest re-anchored with it" = c("aft", "aft.var"),
  "hetero|setWeights|gaussian" = "hetero",
  "hetero|setWeights|aft" = "aft.var",
  "hetero|getLatents|aft" = "aft.var"
)

## ---- 7. The checks ---------------------------------------------------------

checks <- list()
# act(setupResult, seed) runs inside attempt(); setup(seed) does not, so a
# fixture that fails to build is a harness failure rather than a refusal -
# unless the claim is S (then the failure is the disagreement) or the check
# says the fixture's construction is itself the claim (setupClaim)
addCheck <- function(
  source,
  label,
  expect,
  act,
  setup = NULL,
  anchor = NULL,
  setupClaim = FALSE
) {
  checks[[length(checks) + 1L]] <<- list(
    source = source,
    label = label,
    expect = expect,
    act = act,
    setup = setup,
    anchor = anchor,
    setupClaim = setupClaim
  )
}
addMutation <- function(source, label, expect, variant, probe, anchor = NULL) {
  force(variant)
  force(probe)
  addCheck(
    source,
    label,
    expect,
    act = function(s, seed) mutate[[probe]](s, variant),
    setup = function(seed) buildBase(variant, seed),
    anchor = anchor
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
    return(invisible()) # no sampler, no concept: the row's footnote [f10] says why
  }
  if (cap %in% c("xbart", "flatc")) {
    if (cap == "xbart") {
      addCheck(
        "table2",
        label,
        cl$token,
        act = function(s, seed) doXbart(fam, seed)
      )
    } else {
      addCheck(
        "table2",
        label,
        cl$token,
        act = function(s, seed) doFlatC(s),
        setup = function(seed) buildBase(fam, seed),
        setupClaim = fam == "hurdle"
      )
    }
    return(invisible())
  }
  if (cl$token == "-") {
    if (cap == "updateScale") {
      # [f7]: ignored rather than refused where the transform is fixed
      addMutation("table1", label, "S", fam, cap)
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
  for (v in variants) {
    addMutation(
      "table1",
      if (length(variants) > 1L) paste0(label, " [", v, "]") else label,
      cl$token,
      v,
      cap
    )
  }
}))

# -- prose: construction reach
reachAnchor <- "gaussian, student, probit, logistic and aft reach all three"
bartBTAnchor <- "are out of `bartBT()`'s reach entirely"
invisible(lapply(
  c("gaussian", "student", "probit", "logistic", "aft"),
  function(fam) {
    addCheck(
      "prose",
      sprintf("reach %-8s bartBT", fam),
      "S",
      function(s, seed) doBartBT(fam, seed),
      anchor = reachAnchor
    )
    addCheck(
      "prose",
      sprintf("reach %-8s bart", fam),
      "S",
      function(s, seed) doBart(fam, seed),
      anchor = reachAnchor
    )
    addMutation(
      "prose",
      sprintf("reach %-8s dbarts", fam),
      "S",
      fam,
      "getLatents",
      anchor = reachAnchor
    )
  }
))
invisible(lapply(
  c("ordinal", "nbinom", "multinom", "hurdle", "hazard", "bcf"),
  function(fam) {
    addCheck(
      "prose",
      sprintf("reach %-8s bartBT", fam),
      "M",
      function(s, seed) doBartBT(fam, seed),
      anchor = bartBTAnchor
    )
  }
))
addCheck(
  "prose",
  "reach hetero   bartBT",
  "M",
  function(s, seed) {
    d <- mkXY(seed)
    bartBT(d$x, d$y, variance = varForest())
  },
  anchor = bartBTAnchor
)
invisible(lapply(
  c(
    "ordinal",
    "nbinom",
    "multinom",
    "hurdle",
    "hazard",
    "bcf",
    "hetero"
  ),
  function(fam) {
    addCheck(
      "prose",
      sprintf("reach %-8s bart", fam),
      "S",
      function(s, seed) doBart(fam, seed),
      anchor = "`bart()` ships all six"
    )
  }
))
invisible(lapply(
  c("ordinal", "nbinom", "multinom", "hazard", "bcf", "hetero"),
  function(fam) {
    addCheck(
      "prose",
      sprintf("reach %-8s dbarts", fam),
      "S",
      function(s, seed) NULL,
      setup = function(seed) buildBase(fam, seed),
      anchor = "`dbarts()` + R5 ships every one but hurdle"
    )
  }
))
addCheck(
  "prose",
  "reach hurdle   dbarts",
  "R",
  function(s, seed) buildBase("hurdle", seed),
  anchor = "`dbarts()` + R5 ships every one but hurdle"
)
addCheck(
  "prose",
  "reach twopart  bart",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    do.call(bart, c(list(d$x, d$y, family = "twopart"), doorArgs()))
  },
  anchor = "is not a family token at either door"
)
addCheck(
  "prose",
  "bartBT factor of three levels",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    bartBT(d$x, factor(d$label), ntree = 3L, ndpost = 1L, nskip = 0L)
  },
  anchor = "a factor of three or more levels"
)

# -- prose: forests = admits gaussian/probit/logistic only; variance forest
# requires gaussian or aft (the "extra:variance" threading that
# docs/design/aft-variance-forest.md cites is the buildBase(extra = list(variance
# = ...)) calls below)
invisible(lapply(c("aft", "ordinal", "nbinom"), function(fam) {
  addCheck(
    "prose",
    sprintf("forests= %-8s dbarts", fam),
    "R",
    function(s, seed) {
      d <- mkXY(seed)
      y <- switch(fam, aft = cbind(d$time, d$status), ordinal = d$yOrd, d$yNb)
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
    anchor = "A `forests =` fit resolves gaussian, probit or logistic only"
  )
}))
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
  setup = NULL,
  anchor = "under an identical gate"
)
invisible(lapply(c("gaussian", "aft"), function(fam) {
  v <- if (fam == "aft") "aft.var" else "hetero"
  addCheck(
    "prose",
    sprintf("variance= %-8s", fam),
    "S",
    function(s, seed) NULL,
    setup = function(seed) buildBase(v, seed),
    anchor = "requires `family = \"gaussian\"` or `\"aft\"`"
  )
}))
invisible(lapply(
  c("probit", "logistic", "ordinal", "nbinom", "student"),
  function(fam) {
    addCheck(
      "prose",
      sprintf("variance= %-8s", fam),
      "R",
      function(s, seed) buildBase(fam, seed, list(variance = varForest())),
      anchor = if (fam == "student") {
        "it refuses Student-t residuals"
      } else {
        "The other four single-forest families"
      }
    )
  }
))
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
  anchor = "not supported with monotone constraints"
)

# -- prose: DART refused for both couplings
dartAnchor <- "DART is refused for both multi-forest couplings by name"
addCheck(
  "prose",
  "dart bcf",
  "R",
  function(s, seed) {
    buildBase("bcf", seed, list(tree.prior = dbarts:::dart()))
  },
  anchor = dartAnchor
)
addCheck(
  "prose",
  "dart multinom",
  "R",
  function(s, seed) {
    d <- mkXY(seed)
    do.call(
      bart,
      c(bartArgs("multinom", d), doorArgs(tree.prior = dbarts:::dart()))
    )
  },
  anchor = dartAnchor
)

# -- prose: warm start and grow-from-root
warmAnchor <- "Warm start and grow-from-root are refused, as `M` rather than `R`, for the four alternate-family `bart` arcs"
invisible(lapply(c("ordinal", "nbinom", "multinom", "hurdle"), function(fam) {
  addCheck(
    "prose",
    sprintf("grow-from-root %-8s bart", fam),
    "M",
    function(s, seed) {
      d <- mkXY(seed)
      do.call(bart, c(bartArgs(fam, d), doorArgs(n.grow.sweeps = 1L)))
    },
    anchor = warmAnchor
  )
}))
invisible(lapply(
  c("gaussian", "student", "probit", "logistic", "aft", "hazard"),
  function(fam) {
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
      function(s, seed) s$sampler$installTrees(s$donor),
      setup = function(seed) {
        list(
          sampler = buildBase(fam, seed),
          donor = buildBase(fam, seed + 1000L)
        )
      },
      anchor = warmAnchor
    )
  }
))
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
  anchor = "and for the multi-forest donor warm start"
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
invisible(lapply(samplerVariants, function(v) {
  for (probe in c("setResponse", "setOffset", "setPredictor")) {
    blocked <- v == "multinom" && probe != "setPredictor"
    addMutation(
      "prose",
      sprintf("%-8s %s", v, probe),
      if (blocked) "R" else "S",
      v,
      probe,
      anchor = openAnchor
    )
  }
}))

# -- prose: zero weights follow the case-weights column except logistic;
# probit and ordinal take a 0/1 weight vector at creation
zeroAnchor <- "A zero weight, which takes a row out of the likelihood"
invisible(lapply(
  c(
    "gaussian",
    "student",
    "probit",
    "ordinal",
    "hetero",
    "bcf",
    "logistic",
    "bcf.logistic"
  ),
  function(v) {
    addCheck(
      "prose",
      sprintf("zero weight at creation %-8s", v),
      if (grepl("logistic", v, fixed = TRUE)) "R" else "S",
      function(s, seed) {
        buildBase(v, seed, list(weights = c(0, rep(1, n - 1L))))
      },
      anchor = zeroAnchor
    )
  }
))

# -- prose: the universal capabilities
invisible(lapply(
  c(samplerVariants, "aft.var", "bcf.probit", "bcf.logistic"),
  function(v) {
    addMutation(
      "prose",
      sprintf("active-rows mask %-12s", v),
      "S",
      v,
      "activeRows",
      anchor = "the mid-chain active-rows mask"
    )
  }
))
invisible(lapply(
  c(
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
  ),
  function(fam) {
    addCheck(
      "prose",
      sprintf("extract loglik %-8s", fam),
      "S",
      function(s, seed) extract(s, type = "loglik"),
      setup = function(seed) doBart(fam, seed),
      anchor = "`extract(type = \"loglik\")` (every row"
    )
  }
))
invisible(lapply(
  c(
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
  ),
  function(v) {
    coupled <- v %in% c("bcf", "multinom")
    addMutation(
      "prose",
      sprintf("named leaf sd %-8s", v),
      if (coupled) "R" else "S",
      v,
      "leafSd",
      anchor = "A named leaf-prior sd"
    )
  }
))

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

## ---- 8. Run every check, diff against the matrix --------------------------

cat(sprintf(
  "composition-matrix: %d checks (%d table 1, %d table 2, %d prose)\n",
  length(checks),
  sum(vapply(checks, function(ck) ck$source == "table1", logical(1L))),
  sum(vapply(checks, function(ck) ck$source == "table2", logical(1L))),
  sum(vapply(checks, function(ck) ck$source == "prose", logical(1L)))
))

confirmations <- 0L
notes <- character()
disagreements <- character()
harnessFailures <- character()

for (ck in checks) {
  seed <- sum(utf8ToInt(ck$label)) %% 100000L + 1L
  setupResult <- NULL
  if (!is.null(ck$setup)) {
    setupResult <- tryCatch(
      ck$setup(seed),
      error = function(e) e
    )
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
  expect <- ck$expect
  if (nzchar(flipText) && grepl(flipText, ck$label, fixed = TRUE)) {
    expect <- switch(expect, S = "R", R = "S", M = "S")
  }
  agree <- switch(
    expect,
    S = result$verdict == "CONSTRUCTS",
    R = result$verdict == "REFUSES",
    M = result$verdict != "CONSTRUCTS"
  )
  line <- sprintf(
    "[%s] %s: matrix %s, got %s%s",
    ck$source,
    ck$label,
    expect,
    result$verdict,
    if (nzchar(result$detail)) sprintf(" (%s)", result$detail) else ""
  )
  if (verbose) {
    cat(if (agree) "ok   " else "FAIL ", line, "\n")
  }
  if (agree) {
    confirmations <- confirmations + 1L
    if (expect == "M" && result$verdict == "ERRORS") {
      notes <- c(notes, line)
    }
  } else {
    disagreements <- c(disagreements, line)
  }
}

cat(sprintf("\nconfirmations: %d of %d\n", confirmations, length(checks)))
cat(sprintf("\nnotes, generic error on an unbuilt cell (%d):\n", length(notes)))
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
