#!/usr/bin/env Rscript

# Executable composition matrix: for every (family x capability) cell that
# docs/design/feature-matrix.md marks S (shipped) or ? (unverified), builds a
# minimal fixture, actually attempts the composition, and classifies the
# outcome as CONSTRUCTS (fits without error), REFUSES (a named validation
# stop()) or ERRORS (anything else - the finding class: a generic R/engine
# failure rather than a deliberate refusal). The recorded verdict is then
# diffed against the matrix's own claim: an S cell that does not CONSTRUCT is
# a disagreement, a ? cell's verdict is a resolution (reported, never a
# failure), and an R cell that CONSTRUCTS - reachable only if some other
# probe happens to also cover it - would be a disagreement too. Lives in
# benchmarks/R/ rather than tools/: it is a local correctness gate against a
# design document, the same role equivalence.R and bench-sampler.R play
# against a baseline, not a build-freshness check.
#
# Every fixture is tiny (n = 40, a handful of trees, 0 burn-in, 1 sample,
# n.threads = 1L throughout - CRAN's core limit fires before any validation
# this runner cares about) and reseeded per cell, so a run is deterministic
# and fast; nothing here is a statistical check.
#
# Usage: Rscript benchmarks/R/composition-matrix.R
# Exit status is nonzero iff at least one disagreement was recorded.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

scriptDir <- dirname(sub(
  "--file=",
  "",
  grep("--file=", commandArgs(), value = TRUE)
))
matrixPath <- file.path(
  scriptDir,
  "..",
  "..",
  "docs",
  "design",
  "feature-matrix.md"
)

## ---- 1. Parse the matrix's S/? cells out of tables 1-4 --------------------

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
  "grouped",
  "hetero"
)
# transcribed once from the doc's table headers (docs/design/feature-matrix.md
# section titles), used only as dispatch keys and a header sanity check - the
# cell VALUES themselves are read out of the file, not transcribed
tableCols <- list(
  bart.R = c("bart", "bart2", "dbarts5", "rbart_vi", "xbart", "flatc"),
  R.R = c(
    "setResponse",
    "setOffset",
    "updateScale",
    "setPredictor",
    "setWeights",
    "setSigma",
    "testSurface"
  ),
  R.R2 = c(
    "zeroWeightSubset",
    "activeRowsMask",
    "getLatents",
    "pointwiseLoglik",
    "calibration"
  ),
  R.R3 = c(
    "varianceForest",
    "groupedRanef",
    "dart",
    "warmStart",
    "growFromRoot"
  )
)
tableMarkers <- c(
  "`bart2()`",
  "`setResponse`",
  "zero-weight row subset",
  "variance forest"
)
names(tableCols) <- tableMarkers

parseCells <- function(path) {
  lines <- readLines(path)
  legend <- c("S", "R", "P", "M", "-", "?")
  cells <- list()
  for (marker in tableMarkers) {
    hdr <- grep(marker, lines, fixed = TRUE)
    hdr <- hdr[startsWith(trimws(lines[hdr]), "| model |")]
    if (length(hdr) != 1L) {
      stop("composition-matrix: no unique table header for '", marker, "'")
    }
    cols <- tableCols[[marker]]
    body <- lines[seq(hdr + 2L, length.out = length(rowKeys))]
    for (i in seq_along(rowKeys)) {
      raw <- strsplit(body[i], "|", fixed = TRUE)[[1]]
      raw <- trimws(raw[-1L]) # a leading "| " splits off a leading "", never a trailing one
      if (length(raw) != length(cols) + 1L) {
        stop(
          "composition-matrix: '",
          rowKeys[i],
          "' under '",
          marker,
          "' has ",
          length(raw),
          " cells, expected ",
          length(cols) + 1L
        )
      }
      key <- gsub("`", "", raw[1L])
      if (!identical(key, rowKeys[i])) {
        stop(
          "composition-matrix: row ",
          i,
          " under '",
          marker,
          "' reads '",
          key,
          "', expected '",
          rowKeys[i],
          "' - row order drifted"
        )
      }
      for (j in seq_along(cols)) {
        token <- sub("^(\\S+).*$", "\\1", raw[j + 1L])
        if (!token %in% legend) {
          stop(
            "composition-matrix: cell (",
            rowKeys[i],
            ", ",
            cols[j],
            ") = '",
            raw[j + 1L],
            "' has no legend token - ",
            "matrix vocabulary needs a human look"
          )
        }
        if (token %in% c("S", "?")) {
          cells[[length(cells) + 1L]] <- list(
            family = rowKeys[i],
            capability = cols[j],
            claim = token
          )
        }
      }
    }
  }
  cells
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

## ---- 3. Per-family base fixture (the R5 dbartsSampler tables 2-4 mutate) --
## Everything but multinom/hurdle: those two host-shell families carry only
## one or two S/? cells total (table1/table3 and table4 respectively) and are
## special-cased directly in their probes instead of forcing a 13th shape into
## this recipe.

buildBase <- function(family, seed, extra = list()) {
  d <- mkXY(seed)
  args <- switch(
    family,
    gaussian = list(
      d$x,
      d$y,
      test = d$x[1:5, , drop = FALSE],
      control = ctl(seed)
    ),
    student = list(
      d$x,
      d$y,
      test = d$x[1:5, , drop = FALSE],
      resid.dist = dbarts:::student(df = 5),
      control = ctl(seed)
    ),
    probit = list(
      d$x,
      d$yBin,
      test = d$x[1:5, , drop = FALSE],
      family = "probit",
      control = ctl(seed)
    ),
    logistic = list(
      d$x,
      d$yBin,
      test = d$x[1:5, , drop = FALSE],
      family = "logistic",
      control = ctl(seed)
    ),
    ordinal = list(
      d$x,
      d$yOrd,
      test = d$x[1:5, , drop = FALSE],
      family = "ordinal",
      control = ctl(seed)
    ),
    nbinom = list(
      d$x,
      d$yNb,
      test = d$x[1:5, , drop = FALSE],
      family = "nbinom",
      dispersion = 5,
      control = ctl(seed)
    ),
    aft = list(
      d$x,
      cbind(d$time, d$status),
      test = d$x[1:5, , drop = FALSE],
      family = "aft",
      control = ctl(seed)
    ),
    hazard = list(
      d$x,
      survival::Surv(d$time, d$status),
      family = "hazard",
      control = ctl(seed)
    ),
    # no 'test =': a treatment forest refuses test predictors outright
    # (docs/design/bcf.md), and none of the cells this fixture backs need one
    bcf = list(
      d$x,
      d$y,
      forests = list(forest(), forest(basis = ~ factor(d$z))),
      control = ctl(seed)
    ),
    grouped = {
      control <- ctl(seed)
      attr(control, "bartcore.groups") <- list(
        indices = as.integer(d$group),
        n.groups = nlevels(factor(d$group)),
        prior = "cauchy",
        rel.scale = sd(d$y),
        n.steps = 1L
      )
      list(d$x, d$y, test = d$x[1:5, , drop = FALSE], control = control)
    },
    hetero = list(
      d$x,
      d$y,
      test = d$x[1:5, , drop = FALSE],
      variance = varianceForest(n.trees = 3L),
      control = ctl(seed)
    ),
    stop("composition-matrix: no base fixture recipe for '", family, "'")
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

## ---- 5. Table 1: construction surfaces ------------------------------------

probeBart1 <- function(family, seed) {
  d <- mkXY(seed)
  a <- function(...) {
    list(
      ntree = 3L,
      ndpost = 1L,
      nskip = 0L,
      nchain = 1L,
      nthread = 1L,
      verbose = FALSE,
      ...
    )
  }
  switch(
    family,
    gaussian = do.call(bart, c(list(d$x, d$y), a())),
    student = do.call(
      bart,
      c(list(d$x, d$y, resid.dist = dbarts:::student(df = 5)), a())
    ),
    probit = do.call(bart, c(list(d$x, d$yBin), a())),
    logistic = do.call(bart, c(list(d$x, d$yBin, family = "logistic"), a())),
    aft = do.call(
      bart,
      c(list(d$x, cbind(d$time, d$status), family = "aft"), a())
    ),
    stop("composition-matrix: no bart() recipe for '", family, "'")
  )
}

# the bart2() argument recipe shared by probeBart2 (table 1) and
# pointwiseLoglik (table 3) for every family that reaches bart2() plainly;
# NULL for the families each caller must still handle on its own (bcf's
# forest() term, ordinal/nbinom/multinom/hurdle for probeBart2 only)
bart2Args <- function(family, d) {
  switch(
    family,
    gaussian = list(d$x, d$y),
    student = list(d$x, d$y, resid.dist = dbarts:::student(df = 5)),
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
    hetero = list(d$x, d$y, variance = varianceForest(n.trees = 3L)),
    NULL
  )
}

probeBart2 <- function(family, seed) {
  d <- mkXY(seed)
  a <- function(...) {
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
  if (family == "bcf") {
    df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L], z = d$z, y = d$y)
    return(do.call(bart2, c(list(y ~ x1 + x2 + z:forest(x1 + x2), df), a())))
  }
  args <- bart2Args(family, d)
  if (is.null(args)) {
    stop("composition-matrix: no bart2() recipe for '", family, "'")
  }
  do.call(bart2, c(args, a()))
}

probeRbart <- function(family, seed) {
  d <- mkXY(seed)
  grp <- factor(d$group)
  a <- function(...) {
    list(
      n.trees = 3L,
      n.burn = 0L,
      n.samples = 1L,
      n.chains = 1L,
      n.threads = 1L,
      n.thin = 1L,
      verbose = FALSE,
      group.by = grp,
      ...
    )
  }
  switch(
    family,
    gaussian = do.call(rbart_vi, c(list(d$x, d$y), a())),
    probit = do.call(rbart_vi, c(list(d$x, d$yBin), a())),
    aft = {
      df <- data.frame(
        x1 = d$x[, 1L],
        x2 = d$x[, 2L],
        time = d$time,
        status = d$status
      )
      do.call(
        rbart_vi,
        c(list(survival::Surv(time, status) ~ x1 + x2, df, family = "aft"), a())
      )
    },
    grouped = do.call(rbart_vi, c(list(d$x, d$y), a())),
    stop("composition-matrix: no rbart_vi() recipe for '", family, "'")
  )
}

probeXbart <- function(family, seed) {
  d <- mkXY(seed)
  df <- data.frame(x1 = d$x[, 1L], x2 = d$x[, 2L])
  a <- function(...) {
    list(
      n.samples = 2L,
      n.reps = 1L,
      n.burn = c(1L, 1L),
      n.trees = 3L,
      n.threads = 1L,
      verbose = FALSE,
      ...
    )
  }
  switch(
    family,
    gaussian = do.call(xbart, c(list(y ~ x1 + x2, cbind(df, y = d$y)), a())),
    probit = do.call(
      xbart,
      c(list(y ~ x1 + x2, cbind(df, y = d$yBin), family = "probit"), a())
    ),
    logistic = do.call(
      xbart,
      c(list(y ~ x1 + x2, cbind(df, y = d$yBin), family = "logistic"), a())
    ),
    stop("composition-matrix: no xbart() recipe for '", family, "'")
  )
}

# the flat C dbarts.h route: the same compiled consumer test-capi.R drives,
# resolved through R_GetCCallable exactly as a LinkingTo package would. Reuses
# whichever base fixture the family already needed for tables 2-4 - the R5
# sampler's $control/$model/$data slots are the flat API's own creation
# arguments, so this is not a second construction path, just a second entry
# point onto the first one.
capiDll <- local({
  consumerSource <- system.file(
    "tinytest",
    "capi",
    "consumer.c",
    package = "dbarts"
  )
  if (!nzchar(consumerSource)) {
    return(NULL)
  }
  buildDir <- tempfile("capi")
  dir.create(buildDir)
  file.copy(consumerSource, file.path(buildDir, "consumer.c"))
  includeDir <- system.file("include", package = "dbarts")
  writeLines(
    sprintf('PKG_CPPFLAGS = -I"%s"', includeDir),
    file.path(buildDir, "Makevars")
  )
  owd <- setwd(buildDir)
  on.exit(setwd(owd))
  system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "consumer.c"),
    stdout = FALSE,
    stderr = FALSE
  )
  lib <- file.path(buildDir, paste0("consumer", .Platform$dynlib.ext))
  if (!file.exists(lib)) NULL else dyn.load(lib)
})

probeFlatC <- function(family, seed) {
  if (is.null(capiDll)) {
    stop("could not compile the flat C API consumer")
  }
  sampler <- buildBase(family, seed)
  # every other family resolves off model@family / control attributes with an
  # empty family string; aft's survival status attribute needs the explicit
  # token or creation reads it against a defaulted, non-aft family
  token <- if (family == "aft") "aft" else ""
  .Call(
    getNativeSymbolInfo("capi_create", capiDll),
    sampler$control,
    sampler$model,
    sampler$data,
    token
  )
}

table1Probes <- list(
  bart = probeBart1,
  bart2 = probeBart2,
  dbarts5 = function(family, seed) buildBase(family, seed),
  rbart_vi = probeRbart,
  xbart = probeXbart,
  flatc = probeFlatC
)

## ---- 6. Tables 2-4: mutation channels, subsetting, composition -----------

# multinom and hurdle carry only one S/? cell apiece outside table 1
# (multinom: table3 active-rows mask; hurdle: table4 DART); special-cased
# rather than forced through the shared R5 recipe above.
multinomActiveRows <- function(seed) {
  d <- mkXY(seed)
  host <- dbarts(d$x, as.double(d$label), control = ctl(seed))
  bc <- dbarts:::bartcoreMultinomialSampler(host, d$label - 1L, K = 3L)
  bartcoreSetActiveRows(bc, c(rep(1, n - 1L), 0))
}
hurdleDart <- function(seed) {
  d <- mkXY(seed)
  bart2(
    d$x,
    ifelse(d$yBin == 1L, exp(d$x[, 1L]) + 0.1, 0),
    family = "hurdle.lognormal",
    dart = TRUE,
    n.trees = 3L,
    n.burn = 0L,
    n.samples = 1L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}

# pointwise loglik (table 3) dispatches on a bart-classed FIT object
# (extract.bart / extract.rbart), not the raw R5 sampler every other probe
# below mutates - its own small per-family fixture builder
pointwiseLoglik <- function(family, seed) {
  d <- mkXY(seed)
  a <- function(...) {
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
  if (family == "grouped") {
    fit <- do.call(
      rbart_vi,
      c(list(d$x, d$y, group.by = factor(d$group), n.thin = 1L), a())
    )
    return(extract(fit, type = "loglik"))
  }
  args <- bart2Args(family, d)
  if (is.null(args)) {
    stop("composition-matrix: no pointwise-loglik fixture for '", family, "'")
  }
  extract(do.call(bart2, c(args, a())), type = "loglik")
}

# generic capability probes shared by every family the base R5 recipe covers.
# Every probe derives its own row count from the sampler rather than closing
# over the fixture's nominal n: the hazard family expands (x, time, status)
# into a longer person-period design at creation (docs/design/survival.md),
# so its sampler's actual shape disagrees with the un-expanded fixture.
mutate <- list(
  setResponse = function(s, d) s$setResponse(s$data@y),
  setOffset = function(s, d) s$setOffset(rep(0, length(s$data@y))),
  updateScale = function(s, d) {
    s$setOffset(rep(0, length(s$data@y)), updateScale = TRUE)
  },
  setPredictor = function(s, d) s$setPredictor(s$data@x),
  setWeights = function(s, d) s$setWeights(runif(length(s$data@y), 0.5, 1.5)),
  setSigma = function(s, d) s$setSigma(1.0),
  testSurface = function(s, d) {
    x <- s$data@x
    s$setTestPredictor(x[seq_len(min(5L, nrow(x))), , drop = FALSE])
  },
  zeroWeightSubset = "extra:weights",
  activeRowsMask = function(s, d) {
    nr <- length(s$data@y)
    s$setActiveRows(c(rep(1, nr - 1L), 0))
  },
  getLatents = function(s, d) s$getLatents(),
  calibration = function(s, d) s$setCalibration(prior.scale = 1.2),
  varianceForest = "extra:variance",
  groupedRanef = "extra:groups",
  dart = "extra:dart",
  warmStart = function(s, d) {
    donor <- buildBase(d$family, d$seed + 1000L)
    s$installTrees(donor)
  },
  growFromRoot = function(s, d) s$growFromRoot(1L)
)

runProbe <- function(family, capability, seed) {
  if (family == "multinom" && capability == "activeRowsMask") {
    return(attempt(multinomActiveRows(seed)))
  }
  if (family == "hurdle" && capability == "dart") {
    return(attempt(hurdleDart(seed)))
  }
  if (capability %in% names(table1Probes)) {
    return(attempt(table1Probes[[capability]](family, seed)))
  }
  if (capability == "pointwiseLoglik") {
    return(attempt(pointwiseLoglik(family, seed)))
  }
  probe <- mutate[[capability]]
  if (identical(probe, "extra:weights")) {
    return(attempt(buildBase(
      family,
      seed,
      list(weights = c(0, rep(1, n - 1L)))
    )))
  }
  if (identical(probe, "extra:variance")) {
    return(attempt(buildBase(
      family,
      seed,
      list(variance = varianceForest(n.trees = 3L))
    )))
  }
  if (identical(probe, "extra:groups")) {
    grp <- mkXY(seed)$group
    control <- ctl(seed)
    attr(control, "bartcore.groups") <- list(
      indices = as.integer(grp),
      n.groups = nlevels(factor(grp)),
      prior = "cauchy",
      rel.scale = 1,
      n.steps = 1L
    )
    return(attempt(buildBase(family, seed, list(control = control))))
  }
  if (identical(probe, "extra:dart")) {
    return(attempt(buildBase(family, seed, list(tree.prior = dbarts:::dart()))))
  }
  d <- mkXY(seed)
  d$family <- family
  d$seed <- seed
  sampler <- buildBase(family, seed)
  attempt(probe(sampler, d))
}

## ---- 7. Run every S/? cell, diff against the matrix's claim --------------

cells <- parseCells(matrixPath)
cat(sprintf(
  "composition-matrix: %d cells claimed S or ? across %d families\n",
  length(cells),
  length(rowKeys)
))

confirmations <- 0L
resolutions <- list()
disagreements <- list()

for (cell in cells) {
  seed <- sum(utf8ToInt(paste(cell$family, cell$capability))) %% 100000L + 1L
  result <- tryCatch(
    runProbe(cell$family, cell$capability, seed),
    error = function(e) {
      list(verdict = "ERRORS", detail = paste("harness:", conditionMessage(e)))
    }
  )
  label <- sprintf("%-10s %-16s", cell$family, cell$capability)
  if (cell$claim == "?") {
    suffix <- if (nzchar(result$detail)) sprintf(" (%s)", result$detail) else ""
    resolutions[[length(resolutions) + 1L]] <- sprintf(
      "%s -> %s%s",
      label,
      result$verdict,
      suffix
    )
  } else if (result$verdict != "CONSTRUCTS") {
    disagreements[[length(disagreements) + 1L]] <- sprintf(
      "%s claimed S, got %s (%s)",
      label,
      result$verdict,
      result$detail
    )
  } else {
    confirmations <- confirmations + 1L
  }
}

cat(sprintf("\nconfirmations: %d\n", confirmations))
cat(sprintf("\n? resolutions (%d):\n", length(resolutions)))
for (line in resolutions) {
  cat("  ", line, "\n")
}
cat(sprintf("\nDISAGREEMENTS (%d):\n", length(disagreements)))
for (line in disagreements) {
  cat("  ", line, "\n")
}

quit(status = if (length(disagreements) > 0L) 1L else 0L)
