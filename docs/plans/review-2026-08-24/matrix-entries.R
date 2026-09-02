#!/usr/bin/env Rscript
# MATRIX REVIEWER 1 extension: the fitting entries x EVERY family token x the
# conditioning columns docs/design/feature-matrix.md + public-surface.md name
# as accepted or refused. Companion to matrix.R (which crossed conditioning
# against ONE representative family per entry); this closes that coarsening.
#
#   R_LIBS=<lib> Rscript docs/plans/review-2026-08-24/matrix-entries.R [outDir]
#
# Emits <outDir>/matrix-entries-grid.csv
#   entry,family,conditioning,outcome,detail

suppressPackageStartupMessages({
  library(dbarts)
  library(survival)
})

cliArgs <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(cliArgs) >= 1) cliArgs[[1]] else "docs/plans/review-2026-08-24"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
csvPath <- file.path(outDir, "matrix-entries-grid.csv")
grDevices::pdf(NULL)
options(useFancyQuotes = FALSE)  # keep every recorded message ASCII

t0 <- proc.time()[["elapsed"]]
phase <- function(l) cat(sprintf("[%6.1fs] %s\n", proc.time()[["elapsed"]] - t0, l), file = stderr())

conn <- file(csvPath, open = "wt")
writeLines("entry,family,conditioning,outcome,detail", conn)
escapeCsv <- function(s) {
  s <- as.character(s); s[is.na(s)] <- ""
  s <- gsub('"', '""', s, fixed = TRUE); s <- gsub("[\r\n]+", " ", s)
  paste0('"', s, '"')
}
nCells <- 0L
recordRow <- function(entry, family, conditioning, outcome, detail) {
  writeLines(paste(escapeCsv(entry), escapeCsv(family), escapeCsv(conditioning),
                   escapeCsv(outcome), escapeCsv(detail), sep = ","), conn)
  nCells <<- nCells + 1L
  invisible(NULL)
}

# Same raw-error vocabulary matrix.R uses, so the two grids classify alike.
rawErrorPatterns <- c(
  "^'arg' should be one of", "unused argument", "subscript out of bounds",
  "could not find function", "is missing, with no default",
  "missing value where TRUE/FALSE needed", "non-numeric argument",
  "invalid subscript", "object '.*' not found", "attempt to apply non-function",
  "[$] operator is invalid for atomic vectors", "invalid 'type'", "NA/NaN/Inf",
  "unable to find an inherited method", "no applicable method",
  "incorrect number of dimensions", "non-conformable", "cannot coerce",
  "undefined columns selected", "replacement has .* rows, data has",
  "argument is not interpretable as logical",
  "comparison .* is possible only for", "invalid class .*object",
  "is not a function", "^argument \".*\" is missing"
)
classifyError <- function(msg) {
  for (p in rawErrorPatterns) if (grepl(p, msg, perl = TRUE)) return("error-without-reason")
  "refused"
}
describeValue <- function(v) {
  cls <- paste(class(v), collapse = "/")
  d <- dim(v)
  paste(cls, if (!is.null(d)) paste(d, collapse = "x") else paste0("len", length(v)))
}
runCell <- function(entry, family, conditioning, expr) {
  out <- withCallingHandlers(
    tryCatch(list(ok = TRUE, val = eval(expr, envir = globalenv()), msg = NA_character_),
             error = function(e) list(ok = FALSE, val = NULL, msg = conditionMessage(e))),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  if (out$ok) recordRow(entry, family, conditioning, "accepted", describeValue(out$val))
  else        recordRow(entry, family, conditioning, classifyError(out$msg), out$msg)
  invisible(out)
}

## ------------------------------------------------------------------ fixture
set.seed(20260824)
n <- 40L; nTest <- 12L
x <- matrix(runif(n * 3L), n, 3L); colnames(x) <- c("x1", "x2", "x3")
xTest <- matrix(runif(nTest * 3L), nTest, 3L); colnames(xTest) <- colnames(x)
group <- factor(sample(letters[1:4], n, TRUE))
zTrt <- rbinom(n, 1L, 0.5)

yGaussian <- as.numeric(2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3))
yBinary <- factor(rbinom(n, 1L, plogis(2 * (x[, 1L] - 0.5))), levels = c(0, 1))
etaM <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 0)
pM <- exp(etaM); pM <- pM / rowSums(pM)
labM <- vapply(seq_len(n), function(i) sample.int(3L, 1L, prob = pM[i, ]), integer(1))
yMultinom <- factor(c("a", "b", "c")[labM], levels = c("a", "b", "c"))
yCountsMatrix <- t(rmultinom(n, size = 10L, prob = c(0.4, 0.35, 0.25)))
colnames(yCountsMatrix) <- c("a", "b", "c")
yOrdinal <- ordered(c("lo", "mid", "hi")[1L + (x[, 1L] > 0.33) + (x[, 1L] > 0.66)],
                    levels = c("lo", "mid", "hi"))
yNbinom <- rnbinom(n, size = 5L, mu = exp(0.4 * (x[, 1L] - 0.5)))
ySurv <- survival::Surv(rexp(n, rate = exp(-0.4 * (x[, 1L] - 0.5))), rbinom(n, 1L, 0.85))
yHurdle <- ifelse(runif(n) < 0.3, 0, rlnorm(n, meanlog = 0.2 * (x[, 1L] - 0.5)))

familyY <- list(
  auto = yGaussian, gaussian = yGaussian, probit = yBinary, logistic = yBinary,
  aft = ySurv, multinomial = yMultinom, ordinal = yOrdinal, nbinom = yNbinom,
  hazard = ySurv, hazard.probit = ySurv, hazard.logistic = ySurv,
  hurdle.lognormal = yHurdle, twopart = yHurdle)
allTokens <- names(familyY)

sc <- function(...) do.call(dbartsControl, utils::modifyList(
  list(n.chains = 1L, n.threads = 1L, n.trees = 4L, n.burn = 4L, n.samples = 6L,
       keepTrees = TRUE, keepTrainingFits = TRUE), list(...)))

wPos <- runif(n, 0.5, 1.5)
wCount <- rep(1L, n)          # the logistic (PG trial-count) legal spelling
oNum <- rnorm(n, 0, 0.1)
subsetHalf <- seq_len(n) <= n / 2

# sparse fixture: a data.frame carrying a sparseFactor column, formula route
dfSparse <- data.frame(y = yGaussian, x1 = x[, 1L], x2 = x[, 2L])
dfSparse$sf <- dbarts::sparseFactor(
  factor(sample(c("a", "b", "c"), n, TRUE), levels = c("a", "b", "c")),
  reference = "a")

mkCall <- function(fname, pos, named) as.call(c(list(as.name(fname)), pos, named))

## ------------------------------------------- per-entry base call construction
# Base named arguments carry NO conditioning column: each column is added one
# at a time below, so an entry x family cell's outcome is attributable.
baseNamed <- list(
  dbarts   = function(tok) list(family = tok, control = quote(sc()), verbose = FALSE),
  bart2    = function(tok) list(family = tok, n.trees = 4L, n.samples = 6L, n.burn = 4L,
                                n.chains = 1L, n.threads = 1L, keepTrees = TRUE,
                                keepSampler = TRUE, verbose = FALSE, seed = 1L),
  bart     = function(tok) list(family = tok, ntree = 4L, ndpost = 6L, nskip = 4L,
                                nchain = 1L, nthread = 1L, keeptrees = TRUE,
                                keepsampler = TRUE, verbose = FALSE),
  rbart_vi = function(tok) list(family = tok, group.by = quote(group), n.trees = 4L,
                                n.samples = 6L, n.burn = 4L, n.chains = 1L,
                                n.threads = 1L, keepTrees = TRUE, verbose = FALSE),
  xbart    = function(tok) list(family = tok, n.samples = 5L, n.reps = 2L, n.test = 5L,
                                n.burn = c(4L, 3L), n.trees = 4L, n.threads = 1L,
                                method = "k-fold", verbose = FALSE))

# The call is built over NAMES, not embedded values: an entry that deparses
# its own match.call() (several do) would otherwise be handed a 40x3 literal.
basePos <- function(entry, tok) {
  y <- familyY[[tok]]
  if (entry == "bart" && is.factor(y)) y <- as.numeric(y) - 1
  assign(".yv", y, envir = globalenv())
  list(as.name("x"), as.name(".yv"))
}

buildCall <- function(entry, tok, extra = list()) {
  if (entry == "dbartsSpec") {
    return(bquote({
      dd <- dbartsData(x, familyY[[.(tok)]])
      do.call(dbartsSpec, utils::modifyList(
        list(data = dd, control = sc(), family = .(tok)), .(extra)),
        quote = TRUE)
    }))
  }
  named <- utils::modifyList(baseNamed[[entry]](tok), extra)
  mkCall(entry, basePos(entry, tok), named)
}

entries <- c("dbarts", "bart2", "bart", "rbart_vi", "xbart")

## ------------------- Part D: every family token on every entry (13 x 6 = 78)
phase("Part D: family x entry, full cross including out-of-vocabulary aliases")
accepted <- list()
for (entry in c(entries, "dbartsSpec")) {
  for (tok in allTokens) {
    res <- runCell(entry, tok, "<base>", buildCall(entry, tok))
    accepted[[paste0(entry, ":", tok)]] <- res$ok
  }
}
# dbartsSpec's aft door: the design says a survival fit reaches it through the
# `survival =` formal, not through a Surv response.
runCell("dbartsSpec", "aft", "survival=list(status=)", quote({
  dd <- dbartsData(x, log(ySurv[, 1L]))
  dbartsSpec(data = dd, control = sc(), family = "aft",
             survival = list(status = as.integer(ySurv[, 2L])))
}))
runCell("dbartsSpec", "multinomial", "counts-built dbartsData", quote({
  dd <- dbartsData(x, counts = yCountsMatrix)
  dbartsSpec(data = dd, control = sc(), family = "multinomial")
}))
runCell("dbartsSpec", "ordinal", "ordered-factor dbartsData", quote({
  dd <- dbartsData(x, yOrdinal)
  dbartsSpec(data = dd, control = sc(), family = "ordinal")
}))

## ------------------------------ Part E: conditioning x every accepted family
phase("Part E: conditioning columns crossed against every ACCEPTED family")

# Per-entry spelling for each conditioning column. NULL = the entry has no
# spelling at all; we still probe the canonical one, to capture the wall.
condSpec <- list(
  "weights=positive"       = list(dbarts = list(weights = quote(wPos)), bart2 = list(weights = quote(wPos)),
                                  bart = list(weights = quote(wPos)), rbart_vi = list(weights = quote(wPos)),
                                  xbart = list(weights = quote(wPos)), dbartsSpec = list(weights = quote(wPos))),
  "weights=unit counts"    = list(dbarts = list(weights = quote(wCount)), bart2 = list(weights = quote(wCount)),
                                  bart = list(weights = quote(wCount)), rbart_vi = list(weights = quote(wCount)),
                                  xbart = list(weights = quote(wCount)), dbartsSpec = list(weights = quote(wCount))),
  "offset=numeric"         = list(dbarts = list(offset = quote(oNum)), bart2 = list(offset = quote(oNum)),
                                  bart = list(offset = quote(oNum)), rbart_vi = list(offset = quote(oNum)),
                                  xbart = list(offset = quote(oNum)), dbartsSpec = list(offset = quote(oNum))),
  "offset.category="       = list(dbarts = list(offset.category = quote(rep(0, 3L))),
                                  bart2 = list(offset.category = quote(rep(0, 3L))),
                                  bart = list(offset.category = quote(rep(0, 3L))),
                                  rbart_vi = list(offset.category = quote(rep(0, 3L))),
                                  xbart = list(offset.category = quote(rep(0, 3L))),
                                  dbartsSpec = list(offset.category = quote(rep(0, 3L)))),
  "test=matrix"            = list(dbarts = list(test = quote(xTest)), bart2 = list(test = quote(xTest)),
                                  bart = list(x.test = quote(xTest)), rbart_vi = list(test = quote(xTest)),
                                  xbart = list(test = quote(xTest)), dbartsSpec = list(test = quote(xTest))),
  "subset=half"            = list(dbarts = list(subset = quote(subsetHalf)), bart2 = list(subset = quote(subsetHalf)),
                                  bart = list(subset = quote(subsetHalf)), rbart_vi = list(subset = quote(subsetHalf)),
                                  xbart = list(subset = quote(subsetHalf)), dbartsSpec = list(subset = quote(subsetHalf))),
  "resid.dist=student()"   = list(dbarts = list(resid.dist = quote(student())), bart2 = list(resid.dist = quote(student())),
                                  bart = list(resid.dist = quote(student())), rbart_vi = list(resid.dist = quote(student())),
                                  xbart = list(resid.dist = quote(student())), dbartsSpec = list(resid.dist = quote(student()))),
  "variance=~x1"           = list(dbarts = list(variance = quote(~x1)), bart2 = list(variance = quote(~x1)),
                                  bart = list(variance = quote(~x1)), rbart_vi = list(variance = quote(~x1)),
                                  xbart = list(variance = quote(~x1)), dbartsSpec = list(variance = quote(~x1))),
  "forests=list(forest(),forest(basis))" =
                             list(dbarts = list(forests = quote(list(forest(), forest(basis = ~zTrt)))),
                                  bart2 = list(forests = quote(list(forest(), forest(basis = ~zTrt)))),
                                  bart = list(forests = quote(list(forest(), forest(basis = ~zTrt)))),
                                  rbart_vi = list(forests = quote(list(forest(), forest(basis = ~zTrt)))),
                                  xbart = list(forests = quote(list(forest(), forest(basis = ~zTrt)))),
                                  dbartsSpec = list(forests = quote(list(forest(), forest(basis = ~zTrt))))),
  "n.chains=2"             = list(dbarts = list(control = quote(sc(n.chains = 2L))), bart2 = list(n.chains = 2L),
                                  bart = list(nchain = 2L), rbart_vi = list(n.chains = 2L),
                                  xbart = list(n.chains = 2L), dbartsSpec = list(control = quote(sc(n.chains = 2L)))),
  "keepTrees=FALSE"        = list(dbarts = list(control = quote(sc(keepTrees = FALSE))), bart2 = list(keepTrees = FALSE),
                                  bart = list(keeptrees = FALSE), rbart_vi = list(keepTrees = FALSE),
                                  xbart = list(keepTrees = FALSE), dbartsSpec = list(control = quote(sc(keepTrees = FALSE)))),
  "samplerOnly=TRUE"       = list(dbarts = list(samplerOnly = TRUE), bart2 = list(samplerOnly = TRUE),
                                  bart = list(sampleronly = TRUE), rbart_vi = list(samplerOnly = TRUE),
                                  xbart = list(samplerOnly = TRUE), dbartsSpec = list(samplerOnly = TRUE)))

for (cond in names(condSpec)) {
  for (entry in c(entries, "dbartsSpec")) {
    extra <- condSpec[[cond]][[entry]]
    for (tok in allTokens) {
      if (!isTRUE(accepted[[paste0(entry, ":", tok)]])) next  # base already refused
      if (nzchar(Sys.getenv("MXR1_TRACE"))) cat("TRACE", entry, tok, cond, "\n", file = stderr())
      runCell(entry, tok, cond, buildCall(entry, tok, extra))
    }
  }
}

## ---- sparse inputs (formula route; the sparseFactor container is the design's
## own sparse door). Gaussian response only - the container is a PREDICTOR axis.
phase("Part E2: sparse predictor container, every entry")
runCell("dbarts", "gaussian", "sparse predictor (sparseFactor col)",
        quote(dbarts(y ~ x1 + x2 + sf, dfSparse, control = sc(), verbose = FALSE)))
runCell("bart2", "gaussian", "sparse predictor (sparseFactor col)",
        quote(bart2(y ~ x1 + x2 + sf, dfSparse, n.trees = 4L, n.samples = 6L,
                    n.burn = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE, seed = 1L)))
runCell("bart", "gaussian", "sparse predictor (sparseFactor col)",
        quote(bart(dfSparse[, c("x1", "x2", "sf")], dfSparse$y, ntree = 4L,
                   ndpost = 6L, nskip = 4L, nchain = 1L, nthread = 1L, verbose = FALSE)))
runCell("rbart_vi", "gaussian", "sparse predictor (sparseFactor col)",
        quote(rbart_vi(y ~ x1 + x2 + sf, dfSparse, group.by = group, n.trees = 4L,
                       n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE)))
runCell("xbart", "gaussian", "sparse predictor (sparseFactor col)",
        quote(xbart(y ~ x1 + x2 + sf, dfSparse, n.samples = 5L, n.reps = 2L,
                    n.test = 5L, n.burn = c(4L, 3L), n.trees = 4L, n.threads = 1L,
                    method = "k-fold", verbose = FALSE)))
runCell("dbartsData", "gaussian", "sparse predictor (sparseFactor col)",
        quote(dbartsData(y ~ x1 + x2 + sf, dfSparse)))

## bart2's forest() term route (the design's ONLY bart2 door to K forests)
phase("Part E3: bart2 forest() term route, every family")
dfF <- data.frame(y = yGaussian, x1 = x[, 1L], x2 = x[, 2L], z = zTrt)
for (tok in allTokens) {
  dfF$y <- switch(tok,
    probit = , logistic = as.numeric(yBinary) - 1,
    ordinal = yOrdinal, nbinom = yNbinom, multinomial = yMultinom,
    hurdle.lognormal = , twopart = yHurdle, yGaussian)
  runCell("bart2", tok, "forest() term in formula",
          bquote(bart2(y ~ x1 + z:forest(x2), dfF, family = .(tok),
                       n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
                       n.threads = 1L, verbose = FALSE, seed = 1L)))
}

## ------------------------------------- Part F: resolvers, ladders, aliases
phase("Part F: resolver / ladder / alias probes")
rp <- function(label, expr) runCell("resolver", "-", label, expr)

# 1. augmentationFamily vs resolveFamily vocabulary (the C++ disagreement)
for (tok in c("gaussian", "probit", "logistic", "ordinal", "aft", "nbinom", "student"))
  rp(paste0("dbartsDrawLatents(family=", tok, ")"), bquote(
    dbartsDrawLatents(family = .(tok), y = as.numeric(yBinary) - 1,
                      fit = rep(0, n), weights = if (.(tok) == "logistic") wCount else NULL,
                      sigma = if (.(tok) %in% c("aft", "student")) 1 else NULL,
                      dispersion = if (.(tok) == "nbinom") 2 else NULL,
                      cutpoints = if (.(tok) == "ordinal") c(0, 1) else NULL,
                      df = if (.(tok) == "student") 3 else NULL)))

# 2. monotone spellings
for (sp in c("increasing", "decreasing", "inc", "dec", "none", "up"))
  rp(paste0("dbarts(monotone = list(x1 = \"", sp, "\"))"), bquote(
    dbarts(x, yGaussian, control = sc(), verbose = FALSE,
           monotone = list(x1 = .(sp)))))
for (sp in c("increasing", "inc"))
  rp(paste0("bart2(monotone = ~ ", sp, "(x1))"), bquote(
    bart2(yGaussian ~ x1 + x2, data.frame(yGaussian = yGaussian, x1 = x[,1], x2 = x[,2]),
          monotone = .(as.call(list(as.name("~"), as.call(list(as.name(sp), as.name("x1")))))),
          n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
          n.threads = 1L, verbose = FALSE, seed = 1L)))

# 3. makeind(all=)
rp("makeind(df, all = TRUE)",  quote(makeind(data.frame(f = factor(c("a","b","a"))), all = TRUE)))
rp("makeind(df, all = FALSE)", quote(makeind(data.frame(f = factor(c("a","b","a"))), all = FALSE)))

# 4. dbartsValidateComposition / samplePriorPredictive across families
for (tok in allTokens)
  rp(paste0("dbartsValidateComposition(family=", tok, ")"), bquote(
    dbartsValidateComposition(family = .(tok))))

# 5. amplitude prior scale default under each K-forest family
for (tok in c("gaussian", "probit", "logistic", "aft", "ordinal", "nbinom", "multinomial"))
  rp(paste0("dbarts(forests=, family=", tok, ")"), bquote(
    dbarts(x, familyY[[.(tok)]], family = .(tok), control = sc(), verbose = FALSE,
           forests = list(forest(), forest(basis = ~zTrt)))))

# 6. bart()'s own-class redirect on every token outside its vocabulary
for (tok in setdiff(allTokens, c("auto", "logistic", "aft")))
  rp(paste0("bart(family=\"", tok, "\") redirect"), bquote(
    bart(x, yGaussian, family = .(tok), ntree = 4L, ndpost = 6L, nskip = 4L,
         nchain = 1L, nthread = 1L, verbose = FALSE)))
for (tok in setdiff(allTokens, c("auto", "gaussian", "aft")))
  rp(paste0("rbart_vi(family=\"", tok, "\") redirect"), bquote(
    rbart_vi(x, yGaussian, group.by = group, family = .(tok), n.trees = 4L,
             n.samples = 6L, n.burn = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE)))
for (tok in setdiff(allTokens, c("auto", "gaussian", "probit", "logistic")))
  rp(paste0("xbart(family=\"", tok, "\") redirect"), bquote(
    xbart(x, yGaussian, family = .(tok), n.samples = 5L, n.reps = 2L, n.test = 5L,
          n.burn = c(4L, 3L), n.trees = 4L, n.threads = 1L, method = "k-fold",
          verbose = FALSE)))
for (tok in setdiff(allTokens, eval(formals(dbartsSpec)$family)))
  rp(paste0("dbartsSpec(family=\"", tok, "\") redirect"), bquote({
    dd <- dbartsData(x, yGaussian)
    dbartsSpec(data = dd, control = sc(), family = .(tok))
  }))

# 7. the entry-pair agreement probes the prior grid flagged
rp("dbarts(n.thin = -1) via control", quote(dbarts(x, yGaussian, control = sc(n.thin = -1L), verbose = FALSE)))
rp("bart2(n.thin = -1)",  quote(bart2(x, yGaussian, n.thin = -1L, n.trees = 4L, n.samples = 6L,
                                      n.burn = 4L, n.chains = 1L, n.threads = 1L, verbose = FALSE)))
rp("rbart_vi(n.thin = -1)", quote(rbart_vi(x, yGaussian, group.by = group, n.thin = -1L,
                                           n.trees = 4L, n.samples = 6L, n.burn = 4L,
                                           n.chains = 1L, n.threads = 1L, verbose = FALSE)))
rp("bart(keepevery = -1)", quote(bart(x, yGaussian, keepevery = -1L, ntree = 4L, ndpost = 6L,
                                      nskip = 4L, nchain = 1L, nthread = 1L, verbose = FALSE)))
rp("xbart(n.thin = -1)", quote(xbart(x, yGaussian, n.thin = -1L, n.samples = 5L, n.reps = 2L,
                                     n.test = 5L, n.burn = c(4L, 3L), n.trees = 4L,
                                     n.threads = 1L, method = "k-fold", verbose = FALSE)))
for (e in c("dbarts", "bart2", "rbart_vi", "xbart"))
  rp(paste0(e, "(weights containing NA)"), bquote(
    do.call(.(as.name(e)), c(list(x, yGaussian), utils::modifyList(
      if (.(e) == "dbarts") list(control = sc(), verbose = FALSE)
      else if (.(e) == "xbart") list(n.samples = 5L, n.reps = 2L, n.test = 5L,
                                     n.burn = c(4L, 3L), n.trees = 4L, n.threads = 1L,
                                     method = "k-fold", verbose = FALSE)
      else c(list(n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
                  n.threads = 1L, verbose = FALSE),
             if (.(e) == "rbart_vi") list(group.by = group) else list()),
      list(weights = { w <- runif(n); w[1L] <- NA; w }))))))

# 8. n.threads length-2 (the memo's "3 failures on length-2 input")
for (e in c("dbarts", "bart2", "rbart_vi", "xbart", "bart"))
  rp(paste0(e, "(n.threads = c(1,2))"), bquote(
    switch(.(e),
      dbarts   = dbarts(x, yGaussian, control = sc(n.threads = c(1L, 2L)), verbose = FALSE),
      bart2    = bart2(x, yGaussian, n.threads = c(1L, 2L), n.trees = 4L, n.samples = 6L,
                       n.burn = 4L, n.chains = 1L, verbose = FALSE),
      rbart_vi = rbart_vi(x, yGaussian, group.by = group, n.threads = c(1L, 2L),
                          n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L, verbose = FALSE),
      xbart    = xbart(x, yGaussian, n.threads = c(1L, 2L), n.samples = 5L, n.reps = 2L,
                       n.test = 5L, n.burn = c(4L, 3L), n.trees = 4L, method = "k-fold", verbose = FALSE),
      bart     = bart(x, yGaussian, nthread = c(1L, 2L), ntree = 4L, ndpost = 6L,
                      nskip = 4L, nchain = 1L, verbose = FALSE))))

# 9. offset.category, the only documented door
rp("dbartsData(offset.category=) then dbarts(data=)", quote({
  dd <- dbartsData(x, counts = yCountsMatrix, offset.category = c(0, 0.1, -0.1))
  dbartsSpec(data = dd, control = sc(), family = "multinomial")
}))
rp("bart2(data = <dbartsData with offset.category>)", quote({
  dd <- dbartsData(x, counts = yCountsMatrix, offset.category = c(0, 0.1, -0.1))
  bart2(dd, family = "multinomial", n.trees = 4L, n.samples = 6L, n.burn = 4L,
        n.chains = 1L, n.threads = 1L, verbose = FALSE, seed = 1L)
}))

close(conn)
phase(sprintf("done: %d cells -> %s", nCells, csvPath))
