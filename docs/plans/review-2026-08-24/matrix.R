#!/usr/bin/env Rscript
# Executable consistency matrix, second whole-branch review (bartcore, tip b102e17c).
# Re-run against any installed lib:
#   R_LIBS=<lib> Rscript docs/plans/review-2026-08-24/matrix.R [outDir]
#
# Emits <outDir>/matrix-grid.csv (entry,conditioning,generic,type,sample,outcome,detail)
# and prints phase timings to stderr. matrix-results.md is written by hand from the CSV,
# not by this script.

suppressPackageStartupMessages({
  library(dbarts)
  library(survival)
  # as_draws_array/as_draws_df register into POSTERIOR's own S3 table
  # (R/hooks.R .onLoad), and only fire if posterior is already loaded when
  # dbarts loads, or loads later via its packageEvent hook. Load it here so
  # the class x generic existence check sees what a real posterior-using
  # consumer would.
  if (requireNamespace("posterior", quietly = TRUE)) library(posterior)
})

cliArgs <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(cliArgs) >= 1) cliArgs[[1]] else "docs/plans/review-2026-08-24"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
csvPath <- file.path(outDir, "matrix-grid.csv")

grDevices::pdf(NULL) # discard all plot() output; never open an interactive device

t0 <- proc.time()[["elapsed"]]
phase <- function(label) cat(sprintf("[%6.1fs] %s\n", proc.time()[["elapsed"]] - t0, label), file = stderr())

## ---------------------------------------------------------------- recording

conn <- file(csvPath, open = "wt")
writeLines("entry,conditioning,generic,type,sample,outcome,detail", conn)

escapeCsv <- function(s) {
  s <- as.character(s)
  s[is.na(s)] <- ""
  s <- gsub('"', '""', s, fixed = TRUE)
  s <- gsub("[\r\n]+", " ", s)
  paste0('"', s, '"')
}

recordRow <- function(entry, conditioning, generic, type, sample, outcome, detail) {
  writeLines(
    paste(
      escapeCsv(entry), escapeCsv(conditioning), escapeCsv(generic),
      escapeCsv(type), escapeCsv(sample), escapeCsv(outcome), escapeCsv(detail),
      sep = ","
    ),
    conn
  )
  flush(conn)
  invisible(NULL)
}

# Raw, un-narrated R/S4 errors: no argument named, no reason given - the
# "error without reason" outcome. Anything else that reaches stop() with an
# argument name and an explanation is "refused".
rawErrorPatterns <- c(
  "^'arg' should be one of",
  "unused argument",
  "subscript out of bounds",
  "could not find function",
  "is missing, with no default",
  "missing value where TRUE/FALSE needed",
  "non-numeric argument",
  "invalid subscript",
  "object '.*' not found",
  "attempt to apply non-function",
  "[$] operator is invalid for atomic vectors",
  "invalid 'type'",
  "NA/NaN/Inf",
  "unable to find an inherited method",
  "no applicable method",
  "incorrect number of dimensions",
  "non-conformable",
  "cannot coerce",
  "undefined columns selected",
  "replacement has .* rows, data has",
  "argument is not interpretable as logical",
  "comparison .* is possible only for",
  "unable to find an inherited method",
  "invalid class .*object",
  "is not a function"
)
classifyError <- function(msg) {
  for (p in rawErrorPatterns) if (grepl(p, msg, perl = TRUE)) return("error-without-reason")
  "refused"
}

describeValue <- function(v) {
  cls <- paste(class(v), collapse = "/")
  d <- dim(v)
  shape <- if (!is.null(d)) paste(d, collapse = "x") else paste0("len", length(v))
  dn <- dimnames(v)
  dnStr <- if (!is.null(dn)) {
    paste0(
      "dimnames=[",
      paste(vapply(dn, function(z) {
        if (is.null(z)) return("NULL")
        hz <- head(z, 3)
        paste0("{", paste(hz, collapse = ","), if (length(z) > 3) ",..." else "", "}")
      }, character(1)), collapse = ","),
      "]"
    )
  } else ""
  extraAttrs <- setdiff(names(attributes(v)), c("dim", "dimnames", "class", "names"))
  attrStr <- if (length(extraAttrs) > 0) paste0("attrs=[", paste(extraAttrs, collapse = ","), "]") else ""
  trimws(paste(cls, shape, dnStr, attrStr))
}

# Evaluate expr; record one grid row; return list(ok=, val=) for reuse by later phases.
runCell <- function(entry, conditioning, generic = NA, type = NA, sample = NA, expr) {
  out <- withCallingHandlers(
    tryCatch(
      list(ok = TRUE, val = eval(expr, envir = globalenv()), msg = NA_character_),
      error = function(e) list(ok = FALSE, val = NULL, msg = conditionMessage(e))
    ),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  if (out$ok) {
    recordRow(entry, conditioning, generic, type, sample, "accepted", describeValue(out$val))
  } else {
    recordRow(entry, conditioning, generic, type, sample, classifyError(out$msg), out$msg)
  }
  invisible(out)
}

fitAndKeep <- function(entry, conditioning, expr) {
  res <- runCell(entry, conditioning, expr = expr)
  if (res$ok) res$val else NULL
}

## ------------------------------------------------------------------ fixture
# n ~ 40, few trees, few samples - small enough that the whole grid runs in
# minutes. One shared covariate matrix; one response per family token, built
# to be minimally-viable for that family (not to have good posteriors).

set.seed(20260824)
n <- 40L
nTest <- 12L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("x1", "x2", "x3")
xTest <- matrix(runif(nTest * 3L), nTest, 3L)
colnames(xTest) <- colnames(x)
group <- factor(sample(letters[1:4], n, TRUE))
groupTest <- factor(sample(letters[1:4], nTest, TRUE))

yGaussian <- as.numeric(2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3))
pBin <- plogis(2 * (x[, 1L] - 0.5))
yBinary <- factor(rbinom(n, 1L, pBin), levels = c(0, 1))
etaM <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 0)
pM <- exp(etaM); pM <- pM / rowSums(pM)
labM <- vapply(seq_len(n), function(i) sample.int(3L, 1L, prob = pM[i, ]), integer(1))
yMultinom <- factor(c("a", "b", "c")[labM], levels = c("a", "b", "c"))
yCountsMatrix <- t(rmultinom(n, size = 10L, prob = c(0.4, 0.35, 0.25)))
colnames(yCountsMatrix) <- c("a", "b", "c")
yOrdinal <- ordered(
  c("lo", "mid", "hi")[1L + (x[, 1L] > 0.33) + (x[, 1L] > 0.66)],
  levels = c("lo", "mid", "hi")
)
yNbinom <- rnbinom(n, size = 5L, mu = exp(0.4 * (x[, 1L] - 0.5)))
survTime <- rexp(n, rate = exp(-0.4 * (x[, 1L] - 0.5)))
survStatus <- rbinom(n, 1L, 0.85)
ySurv <- survival::Surv(survTime, survStatus)
yHurdle <- ifelse(runif(n) < 0.3, 0, rlnorm(n, meanlog = 0.2 * (x[, 1L] - 0.5)))

familyY <- list(
  auto = yGaussian, gaussian = yGaussian, probit = yBinary, logistic = yBinary,
  aft = ySurv, multinomial = yMultinom, ordinal = yOrdinal, nbinom = yNbinom,
  hazard = ySurv, hazard.probit = ySurv, hazard.logistic = ySurv,
  hurdle.lognormal = yHurdle, twopart = yHurdle
)

sc <- function(...) {
  defaults <- list(
    n.chains = 1L, n.threads = 1L, n.trees = 4L, n.burn = 4L, n.samples = 6L,
    keepTrees = TRUE, keepTrainingFits = TRUE
  )
  do.call(dbartsControl, utils::modifyList(defaults, list(...)))
}

phase("fixture built; starting Part A (family-token acceptance)")

## ------------------------------------------------ Part A: family acceptance
# Every declared family token for every fitting entry, mechanically pulled
# from the entry's own formal default (so this stays correct if the branch's
# vocabulary changes), plus one deliberately-invalid token per entry.

declaredTokens <- function(fn) eval(formals(fn)$family)

fittedObjects <- new.env() # entry:token -> fitted/constructed object, for Part C

tryToken <- function(entry, token, expr) {
  key <- paste0(entry, ":", token)
  val <- fitAndKeep(entry, paste0("family=", token), expr)
  if (!is.null(val)) assign(key, val, envir = fittedObjects)
  val
}

## dbarts()
for (tok in c(declaredTokens(dbarts), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  tryToken("dbarts", tok, quote({
    s <- dbarts(x, y, family = tok, control = sc(), verbose = FALSE)
    r <- s$run(0L, 6L)
    list(sampler = s, run = r)
  }))
}

## bart2()
for (tok in c(declaredTokens(bart2), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  tryToken("bart2", tok, quote(bart2(
    x, y, family = tok, n.trees = 4L, n.samples = 6L, n.burn = 4L,
    n.chains = 1L, n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE,
    keepSampler = TRUE, test = xTest, verbose = FALSE, seed = 1L
  )))
}

## bart() (BayesTree-compat shim: x.train/y.train, own token vocabulary)
for (tok in c(declaredTokens(bart), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  yv <- if (is.factor(y)) as.numeric(y) - 1 else y
  tryToken("bart", tok, quote(bart(
    x.train = x, y.train = yv, x.test = xTest, ntree = 4L, ndpost = 6L,
    nskip = 4L, nchain = 1L, nthread = 1L, keeptrees = TRUE,
    keeptrainfits = TRUE, keepsampler = TRUE, verbose = FALSE, family = tok
  )))
}

## rbart_vi() (needs a grouping factor)
for (tok in c(declaredTokens(rbart_vi), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  tryToken("rbart_vi", tok, quote(rbart_vi(
    x, y, group.by = group, group.by.test = groupTest, test = xTest,
    n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L,
    keepTrees = TRUE, verbose = FALSE, family = tok
  )))
}

## xbart() (cross-validation; small everything to keep it fast)
for (tok in c(declaredTokens(xbart), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  tryToken("xbart", tok, quote(xbart(
    x, y, n.samples = 5L, n.reps = 2L, n.test = 5L, n.burn = c(4L, 3L),
    n.trees = 4L, n.threads = 1L, method = "k-fold", verbose = FALSE, family = tok
  )))
}

## dbartsSpec() (constructor only; separate token vocabulary)
for (tok in c(declaredTokens(dbartsSpec), "zzz-invalid")) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  key <- paste0("dbartsSpec:", tok)
  dd <- tryCatch(dbartsData(x, y), error = function(e) NULL)
  if (is.null(dd)) {
    recordRow("dbartsSpec", paste0("family=", tok), NA, NA, NA, "error-without-reason", "fixture dbartsData() build failed")
    next
  }
  val <- fitAndKeep("dbartsSpec", paste0("family=", tok), quote(
    dbartsSpec(data = dd, control = sc(), family = tok)
  ))
  if (!is.null(val)) assign(key, val, envir = fittedObjects)
}
# Supplementary: multinomial only reaches a live model through dbartsData's
# 'counts' argument (see data.R), never through a plain factor response -
# confirm the token is reachable at all, on the same data dbarts()/bart2()
# accept directly as a factor via their own (x, y) interface.
ddCounts <- tryCatch(dbartsData(x, counts = yCountsMatrix), error = function(e) NULL)
if (!is.null(ddCounts)) {
  val <- fitAndKeep("dbartsSpec", "family=multinomial (counts-based dbartsData)", quote(
    dbartsSpec(data = ddCounts, control = sc(), family = "multinomial")
  ))
  if (!is.null(val)) assign("dbartsSpec:multinomial.counts", val, envir = fittedObjects)
}

## dbartsData() / dbartsControl(): no family formal at all - record that
## mechanically (not by assertion) by actually passing family= and reading
## the resulting error, once each.
runCell("dbartsData", "family=gaussian (dbartsData has no family formal)",
        expr = quote(dbartsData(x, yGaussian, family = "gaussian")))
runCell("dbartsControl", "family=gaussian (dbartsControl has no family formal)",
        expr = quote(dbartsControl(family = "gaussian")))

phase("Part A done; starting Part B (conditioning columns)")

## ------------------------------------------------- Part B: conditioning axis
# One representative accepted family per entry, crossed with each
# conditioning column one at a time (not a full factorial against Part A -
# see report for what this coarsens away). Each column is tried on every
# entry even when the entry's formals plainly lack it, so the "unused
# argument" wall is captured by execution, not by reading formals.

reprY <- yGaussian

# Every *Col wrapper below merges a per-entry base named-argument list with
# the probe's `extra` via modifyList, so a probe that re-specifies a base key
# (control=, n.chains=, ...) overrides rather than duplicates it; positional
# data args (x, reprY, ...) sit outside the merge and are never touched.

## dbarts
dbCol <- function(label, extra) runCell("dbarts", label, expr = bquote(
  do.call(dbarts, c(list(x, reprY), utils::modifyList(
    list(family = "gaussian", control = sc(), verbose = FALSE), .(extra)
  )))
))
# resid.dist is NSE (like variance/k): must appear as literal unevaluated
# syntax in the actual call (match.call() captures it, re-evaluated inside
# the package's own namespace where student()/gaussian() live unexported) -
# a do.call()-with-pre-evaluated-list probe would evaluate student() in this
# script's frame instead, where it is not visible, so this one call is direct.
runCell("dbarts", "resid.dist=student()", expr = quote(
  dbarts(x, reprY, family = "gaussian", control = sc(), resid.dist = student(), verbose = FALSE)
))
dbCol("variance=~x1", quote(list(variance = ~x1)))
dbCol("weights=valid", quote(list(weights = runif(n, 0.5, 1.5))))
dbCol("weights=NA", quote(list(weights = { w <- runif(n); w[1] <- NA; w })))
dbCol("offset=valid", quote(list(offset = rnorm(n, 0, 0.1))))
dbCol("test=matrix", quote(list(test = xTest)))
dbCol("subset=half", quote(list(subset = seq_len(n) <= n / 2)))
dbCol("control=sc(keepTrees=FALSE)", quote(list(control = sc(keepTrees = FALSE))))
dbCol("samplerOnly=TRUE (no such formal)", quote(list(samplerOnly = TRUE)))
dbCol("n.chains=2-via-formal (no such formal)", quote(list(n.chains = 2L)))
dbCol("combineChains=FALSE (no such formal)", quote(list(combineChains = FALSE)))
dbCol("keepSampler=TRUE (no such formal)", quote(list(keepSampler = TRUE)))
runCell("dbarts", "forests=list(forest())", expr = quote(
  dbarts(x, reprY, family = "gaussian", control = sc(),
         forests = list(forest(n.trees = 3L), forest(basis = x[, 2, drop = FALSE], n.trees = 3L)),
         verbose = FALSE)
))
runCell("dbarts", "family=multinomial + offset.category via dbartsData(counts=)", expr = quote({
  ddm <- dbartsData(x, counts = yCountsMatrix, offset.category = matrix(0, n, 3L))
  dbarts(ddm, family = "multinomial", control = sc(), verbose = FALSE)
}))
runCell("dbarts", "sparse predictor (sparseFactor column)", expr = quote({
  fr <- data.frame(x1 = x[, 1])
  fr$f <- sparseFactor(sample(1:5, n, TRUE), levels = 1:5, reference = 1L)
  dbarts(fr, reprY, control = sc(), verbose = FALSE)
}))
runCell("dbarts", "mixed-matrix predictor", expr = quote({
  mm <- makeModelMatrixFromDataFrame(data.frame(x1 = x[, 1], f = factor(sample(letters[1:3], n, TRUE))))
  dbarts(mm, reprY, control = sc(), verbose = FALSE)
}))

## bart2
b2Col <- function(label, extra) runCell("bart2", label, expr = bquote(
  do.call(bart2, c(list(x, reprY), utils::modifyList(
    list(family = "gaussian", n.trees = 4L, n.samples = 6L, n.burn = 4L,
         n.chains = 1L, n.threads = 1L, verbose = FALSE), .(extra)
  )))
))
runCell("bart2", "resid.dist=student()", expr = quote(
  bart2(x, reprY, family = "gaussian", n.trees = 4L, n.samples = 6L, n.burn = 4L,
        n.chains = 1L, n.threads = 1L, resid.dist = student(), verbose = FALSE)
))
b2Col("variance=~x1", quote(list(variance = ~x1)))
b2Col("weights=valid", quote(list(weights = runif(n, 0.5, 1.5))))
b2Col("weights=NA", quote(list(weights = { w <- runif(n); w[1] <- NA; w })))
b2Col("offset=valid", quote(list(offset = rnorm(n, 0, 0.1))))
b2Col("test=matrix", quote(list(test = xTest)))
b2Col("subset=half", quote(list(subset = seq_len(n) <= n / 2)))
b2Col("keepTrees=FALSE", quote(list(keepTrees = FALSE)))
b2Col("keepSampler=TRUE", quote(list(keepSampler = TRUE)))
b2Col("n.chains=2", quote(list(n.chains = 2L)))
b2Col("combineChains=FALSE", quote(list(n.chains = 2L, combineChains = FALSE)))
b2Col("samplerOnly=TRUE", quote(list(samplerOnly = TRUE)))
runCell("bart2", "forests=list(forest())", expr = quote(
  bart2(x, reprY, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L,
        forests = list(forest(n.trees = 3L), forest(basis = x[, 2, drop = FALSE], n.trees = 3L)),
        verbose = FALSE)
))
runCell("bart2", "offset.category (no such formal, via ...)", expr = quote(
  bart2(x, yMultinom, family = "multinomial", n.trees = 4L, n.samples = 6L, n.burn = 4L,
        offset.category = matrix(0, n, 3L), verbose = FALSE)
))
runCell("bart2", "sparse predictor (sparseFactor column)", expr = quote({
  fr <- data.frame(x1 = x[, 1])
  fr$f <- sparseFactor(sample(1:5, n, TRUE), levels = 1:5, reference = 1L)
  bart2(fr, reprY, n.trees = 4L, n.samples = 6L, n.burn = 4L, verbose = FALSE)
}))

## bart (BayesTree-compat)
bCol <- function(label, extra) runCell("bart", label, expr = bquote(
  do.call(bart, utils::modifyList(
    list(x.train = x, y.train = reprY, x.test = xTest, ntree = 4L, ndpost = 6L,
         nskip = 4L, nchain = 1L, nthread = 1L, verbose = FALSE), .(extra)
  ))
))
runCell("bart", "resid.dist=student()", expr = quote(
  bart(x.train = x, y.train = reprY, x.test = xTest, ntree = 4L, ndpost = 6L, nskip = 4L,
       nchain = 1L, nthread = 1L, resid.dist = student(), verbose = FALSE)
))
bCol("variance=~x1 (no such formal)", quote(list(variance = ~x1)))
bCol("weights=valid", quote(list(weights = runif(n, 0.5, 1.5))))
bCol("weights=NA", quote(list(weights = { w <- runif(n); w[1] <- NA; w })))
bCol("offset (no such formal; binaryOffset exists instead)", quote(list(offset = rnorm(n, 0, 0.1))))
bCol("subset=half", quote(list(subset = seq_len(n) <= n / 2)))
bCol("keeptrees=FALSE", quote(list(keeptrees = FALSE)))
bCol("keepsampler=TRUE", quote(list(keepsampler = TRUE)))
bCol("nchain=2", quote(list(nchain = 2L)))
bCol("combinechains=FALSE", quote(list(nchain = 2L, combinechains = FALSE)))
bCol("sampleronly=TRUE", quote(list(sampleronly = TRUE)))
bCol("forests= (no such formal)", quote(list(forests = list(forest()))))
bCol("n.chains=2 (wrong-case formal name)", quote(list(n.chains = 2L)))

## rbart_vi
rCol <- function(label, extra) runCell("rbart_vi", label, expr = bquote(
  do.call(rbart_vi, c(list(x, reprY), utils::modifyList(
    list(group.by = group, group.by.test = groupTest, test = xTest, n.trees = 4L,
         n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE), .(extra)
  )))
))
runCell("rbart_vi", "resid.dist=student() (no such formal)", expr = quote(
  rbart_vi(x, reprY, group.by = group, group.by.test = groupTest, test = xTest, n.trees = 4L,
           n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L, resid.dist = student(),
           verbose = FALSE)
))
rCol("variance=~x1 (no such formal)", quote(list(variance = ~x1)))
rCol("weights=valid", quote(list(weights = runif(n, 0.5, 1.5))))
rCol("weights=NA", quote(list(weights = { w <- runif(n); w[1] <- NA; w })))
rCol("offset=valid", quote(list(offset = rnorm(n, 0, 0.1))))
rCol("subset=half", quote(list(subset = seq_len(n) <= n / 2)))
rCol("keepTrees=FALSE", quote(list(keepTrees = FALSE)))
rCol("keepSampler=TRUE", quote(list(keepSampler = TRUE)))
rCol("n.chains=2", quote(list(n.chains = 2L)))
rCol("combineChains=FALSE", quote(list(n.chains = 2L, combineChains = FALSE)))
rCol("samplerOnly=TRUE (no such formal)", quote(list(samplerOnly = TRUE)))
rCol("forests= (no such formal)", quote(list(forests = list(forest()))))

## xbart
xCol <- function(label, extra) runCell("xbart", label, expr = bquote(
  do.call(xbart, c(list(x, reprY), utils::modifyList(
    list(n.samples = 5L, n.reps = 2L, n.test = 5L, n.burn = c(4L, 3L),
         n.trees = 4L, n.threads = 1L, method = "k-fold", verbose = FALSE), .(extra)
  )))
))
runCell("xbart", "resid.dist=student() (no such formal)", expr = quote(
  xbart(x, reprY, n.samples = 5L, n.reps = 2L, n.test = 5L, n.burn = c(4L, 3L), n.trees = 4L,
        n.threads = 1L, method = "k-fold", resid.dist = student(), verbose = FALSE)
))
xCol("variance=~x1 (no such formal)", quote(list(variance = ~x1)))
xCol("weights=valid", quote(list(weights = runif(n, 0.5, 1.5))))
xCol("weights=NA", quote(list(weights = { w <- runif(n); w[1] <- NA; w })))
xCol("offset=valid", quote(list(offset = rnorm(n, 0, 0.1))))
xCol("subset=half", quote(list(subset = seq_len(n) <= n / 2)))
xCol("test= (no such formal)", quote(list(test = xTest)))
xCol("keepTrees= (no such formal)", quote(list(keepTrees = TRUE)))
xCol("n.chains= (no such formal)", quote(list(n.chains = 2L)))
xCol("forests= (no such formal)", quote(list(forests = list(forest()))))

phase("Part B done; starting Part C (class x generic existence + execution)")

## ------------------------------------- Part C1: class x generic, mechanical
# Existence is a static fact - check it without needing an object at all.

classesOfInterest <- c(
  "bart", "rbart", "bartMultinomial", "bartOrdinal", "bartNegbin", "bartHurdle",
  "pdbart", "pd2bart", "dbartsSampler"
)
coreGenerics <- c("predict", "extract", "fitted", "residuals", "summary", "print", "plot")
extraGenerics <- c("plotTree", "survivalProbabilities", "as_draws_array", "as_draws_df")

existsMethod <- function(generic, cls) {
  m <- tryCatch(getS3method(generic, cls, optional = TRUE), error = function(e) NULL)
  !is.null(m)
}

for (g in coreGenerics) for (cl in classesOfInterest) {
  ex <- existsMethod(g, cl)
  recordRow("class-x-generic", cl, g, NA, NA, if (ex) "exists" else "no-method", "")
}
for (g in extraGenerics) for (cl in classesOfInterest) {
  ex <- existsMethod(g, cl)
  recordRow("class-x-generic-extra", cl, g, NA, NA, if (ex) "exists" else "no-method", "")
}

phase("Part C1 (existence table) done; building representative fit objects")

## --------------------------------------- Part C2: representative fit objects
# One classed object per (family variant) so every discovered method can
# actually be executed, not just found. Reuses Part A's bart2 fits where the
# class matches, adds a few extra variants Part A didn't need (heteroscedastic
# gaussian, auto-dispatch, pdbart/pd2bart, plain dbarts()+run()).

reps <- new.env()
put <- function(name, val) if (!is.null(val)) assign(name, val, envir = reps)

put("bart.gaussian", get0("bart2:gaussian", envir = fittedObjects))
put("bart.probit", get0("bart2:probit", envir = fittedObjects))
put("bart.logistic", get0("bart2:logistic", envir = fittedObjects))
put("bart.aft", get0("bart2:aft", envir = fittedObjects))
# Part A's bart2 loop passed test= unconditionally, which the hazard/hurdle
# branches refuse outright (discrete-time hazard has no 'test' concept) - so
# those two never landed in fittedObjects. Build hazard reps here without
# test=, matching what Part A already established as the real acceptance path.
put("bart.hazard", fitAndKeep("bart2", "family=hazard (rep, no test=)", quote(
  bart2(x, ySurv, family = "hazard", n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
        n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE, verbose = FALSE, seed = 1L)
)))
put("bart.hazard.logistic", fitAndKeep("bart2", "family=hazard.logistic (rep, no test=)", quote(
  bart2(x, ySurv, family = "hazard.logistic", n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
        n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE, verbose = FALSE, seed = 1L)
)))
put("bartMultinomial", get0("bart2:multinomial", envir = fittedObjects))
put("bartOrdinal", get0("bart2:ordinal", envir = fittedObjects))
put("bartNegbin", get0("bart2:nbinom", envir = fittedObjects))
put("bartHurdle", get0("bart2:hurdle.lognormal", envir = fittedObjects))
put("rbart", get0("rbart_vi:gaussian", envir = fittedObjects))

put("bart.heteroscedastic", fitAndKeep("bart2", "family=gaussian;variance=~x1 (heteroscedastic rep)", quote(
  bart2(x, yGaussian, variance = ~x1, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
        n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE, test = xTest, verbose = FALSE, seed = 1L)
)))
put("bart.student", fitAndKeep("bart2", "family=gaussian;resid.dist=student() (rep)", quote(
  bart2(x, yGaussian, resid.dist = student(), n.trees = 4L, n.samples = 6L, n.burn = 4L,
        n.chains = 1L, n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE, test = xTest,
        verbose = FALSE, seed = 1L)
)))
put("bart.auto.binary", fitAndKeep("bart2", "family=auto (binary y, rep)", quote(
  bart2(x, yBinary, family = "auto", n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
        n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE, test = xTest, verbose = FALSE, seed = 1L)
)))
put("bart.nchains2", fitAndKeep("bart2", "family=gaussian;n.chains=2 (rep)", quote(
  bart2(x, yGaussian, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 2L, n.threads = 1L,
        keepTrees = TRUE, keepTrainingFits = TRUE, test = xTest, verbose = FALSE, seed = 1L)
)))

put("pdbart", fitAndKeep("pdbart", "family=gaussian (rep)", quote(
  pdbart(x, yGaussian, xind = 1:2, ntree = 4L, ndpost = 6L, nskip = 4L, verbose = FALSE)
)))
put("pd2bart", fitAndKeep("pd2bart", "family=gaussian (rep)", quote(
  pd2bart(x, yGaussian, xind1 = 1L, xind2 = 2L, ntree = 4L, ndpost = 6L, nskip = 4L, verbose = FALSE)
)))

put("dbartsSampler.gaussian", {
  wrapped <- get0("dbarts:gaussian", envir = fittedObjects)
  if (is.null(wrapped)) NULL else wrapped$sampler
})
put("dbartsSampler.multinomial", {
  wrapped <- get0("dbarts:multinomial", envir = fittedObjects)
  if (is.null(wrapped)) NULL else wrapped$sampler
})

phase(sprintf("Part C2 done: %d representative objects; starting Part C3 (execution grid)", length(ls(reps))))

## ---------------------------------------------------- Part C3: execution grid
# For every representative object, every existing generic for its class,
# every type= value that generic declares (mechanically read off its own
# formal default - so this tracks the branch, not a hand-typed list), plus
# one deliberately-invalid type probe, times sample = train/test where the
# generic takes a sample argument.

repClass <- function(obj) {
  cl <- class(obj)
  cl[cl %in% classesOfInterest][1]
}

typeValuesFor <- function(fn) {
  f <- formals(fn)
  if (is.null(f) || !("type" %in% names(f))) return(NA_character_)
  v <- tryCatch(eval(f$type), error = function(e) NA_character_)
  if (!is.character(v)) return(NA_character_)
  v
}
sampleValuesFor <- function(fn) {
  f <- formals(fn)
  if (is.null(f) || !("sample" %in% names(f))) return(NA_character_)
  v <- tryCatch(eval(f$sample), error = function(e) NA_character_)
  if (!is.character(v)) return(NA_character_)
  v
}

runGeneric <- function(repName, obj, generic) {
  cls <- repClass(obj)
  if (is.na(cls)) return(invisible(NULL))
  m <- tryCatch(getS3method(generic, cls, optional = TRUE), error = function(e) NULL)
  if (is.null(m)) return(invisible(NULL))

  types <- typeValuesFor(m)
  types <- if (identical(types, NA_character_)) NA_character_ else c(types, "zzz-invalid-type")
  samples <- sampleValuesFor(m)
  samples <- if (identical(samples, NA_character_)) NA_character_ else c(samples, "zzz-invalid-sample")

  fNames <- names(formals(m))
  firstArgName <- fNames[1] # some methods take `object`, print/plot take `x`
  newdataArg <- if ("newdata" %in% fNames) list(newdata = xTest) else list()
  groupByArg <- if ("group.by" %in% fNames) list(group.by = group) else list()
  # plotTree.dbartsSampler(object, ...) forwards ... to a required treeNum
  # the S3 formals don't expose; survivalProbabilities always needs times=.
  # Neither is discoverable by formals() alone, so these are keyed by generic.
  treeNumArg <- if (identical(generic, "plotTree")) list(treeNum = 1L) else list()
  if (identical(generic, "survivalProbabilities")) newdataArg <- c(newdataArg, list(times = c(0.5, 1, 1.5)))

  for (ty in types) for (sm in samples) {
    callArgs <- c(stats::setNames(list(obj), firstArgName), newdataArg, groupByArg, treeNumArg)
    if (!identical(ty, NA_character_)) callArgs$type <- ty
    if (!identical(sm, NA_character_)) callArgs$sample <- sm
    runCell(
      paste0("fit:", repName), paste0("class=", cls),
      generic = generic, type = ty, sample = sm,
      expr = bquote(do.call(.(m), .(callArgs)))
    )
  }
  invisible(NULL)
}

allGenerics <- c(coreGenerics, extraGenerics)
repNames <- ls(reps)
for (rn in repNames) {
  obj <- get(rn, envir = reps)
  for (g in allGenerics) runGeneric(rn, obj, g)
}

phase("Part C3 done; starting Part D (entry-pair alignment probes)")

## ------------------------------------------- Part D: entry-pair alignment
# Same input, multiple entries: dbartsSpec vs dbarts() family resolution
# (man/dbartsSpec.Rd claims identical resolution), and a couple of the
# formal-name/behavior splits Part B already exposed, re-run side by side
# for a direct diff.

specTokens <- declaredTokens(dbartsSpec)
for (tok in specTokens) {
  y <- familyY[[tok]]
  if (is.null(y)) y <- yGaussian
  dd <- tryCatch(dbartsData(x, y), error = function(e) NULL)
  if (is.null(dd)) next
  specRes <- tryCatch(dbartsSpec(data = dd, control = sc(), family = tok), error = function(e) e)
  dbartsRes <- tryCatch(dbarts(x, y, family = tok, control = sc(), verbose = FALSE), error = function(e) e)
  specOk <- !inherits(specRes, "error")
  dbartsOk <- !inherits(dbartsRes, "error")
  specFamily <- if (specOk) tryCatch(specRes@family, error = function(e) tryCatch(specRes$family, error = function(e2) NA)) else conditionMessage(specRes)
  dbartsFamily <- if (dbartsOk) tryCatch(dbartsRes$model@family, error = function(e) NA) else conditionMessage(dbartsRes)
  agree <- specOk == dbartsOk && (!specOk || identical(as.character(specFamily), as.character(dbartsFamily)))
  recordRow(
    "entry-pair:dbartsSpec-vs-dbarts", paste0("family=", tok), "family-resolution", NA, NA,
    if (agree) "agree" else "DISAGREE",
    sprintf("dbartsSpec: ok=%s family=%s | dbarts: ok=%s family=%s", specOk, specFamily, dbartsOk, dbartsFamily)
  )
}

## weights=NA: does dbarts() and bart2() (same underlying validator?) refuse identically?
wNA <- { w <- runif(n); w[1] <- NA; w }
r1 <- tryCatch(dbarts(x, yGaussian, weights = wNA, control = sc(), verbose = FALSE), error = function(e) e)
r2 <- tryCatch(bart2(x, yGaussian, weights = wNA, n.trees = 4L, n.samples = 6L, n.burn = 4L, verbose = FALSE), error = function(e) e)
recordRow(
  "entry-pair:dbarts-vs-bart2", "weights=NA", "weights-validation", NA, NA,
  if (inherits(r1, "error") == inherits(r2, "error") &&
      (!inherits(r1, "error") || identical(conditionMessage(r1), conditionMessage(r2)))) "agree" else "DISAGREE",
  sprintf(
    "dbarts: %s | bart2: %s",
    if (inherits(r1, "error")) conditionMessage(r1) else "accepted",
    if (inherits(r2, "error")) conditionMessage(r2) else "accepted"
  )
)

## thinned-count-negative probe across the three entries that expose n.thin/keepevery
r3 <- tryCatch(bart2(x, yGaussian, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.thin = -1L, verbose = FALSE), error = function(e) e)
r4 <- tryCatch(rbart_vi(x, yGaussian, group.by = group, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.thin = -1L, verbose = FALSE), error = function(e) e)
r5 <- tryCatch(bart(x.train = x, y.train = yGaussian, ntree = 4L, ndpost = 6L, nskip = 4L, keepevery = -1L, verbose = FALSE), error = function(e) e)
recordRow(
  "entry-pair:bart2-vs-rbart_vi-vs-bart", "n.thin/keepevery = -1", "negative-thin-validation", NA, NA,
  if (inherits(r3, "error") && inherits(r4, "error") && inherits(r5, "error")) "agree(all refuse)" else "DISAGREE",
  sprintf(
    "bart2: %s | rbart_vi: %s | bart: %s",
    if (inherits(r3, "error")) conditionMessage(r3) else "ACCEPTED",
    if (inherits(r4, "error")) conditionMessage(r4) else "ACCEPTED",
    if (inherits(r5, "error")) conditionMessage(r5) else "ACCEPTED"
  )
)

phase("Part D done; closing")

close(conn)
phase(sprintf("all done: wrote %s", csvPath))
