#!/usr/bin/env Rscript
# Supplement to matrix.R: axis-1 explicitly names "the dbartsSampler
# reference class (its ~46 own methods; census in
# inst/tinytest/test-host-shell-pins.R)" as one of the fitting entries.
# matrix.R's execution grid only reaches dbartsSampler through the two
# registered S3 generics (extract, plotTree); it never calls the other 39
# substantive RC methods directly. This appends one row per method to the
# same matrix-grid.csv (entry = "dbartsSampler-method", conditioning =
# "gaussian sampler", generic = method name).
#
# Usage: R_LIBS=<lib> Rscript matrix-sampler-census.R [outDir]

suppressPackageStartupMessages({
  library(dbarts)
})

cliArgs <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(cliArgs) >= 1) cliArgs[[1]] else "docs/plans/review-2026-08-24"
csvPath <- file.path(outDir, "matrix-grid.csv")
stopifnot(file.exists(csvPath))

escapeCsv <- function(s) {
  s <- as.character(s)
  s[is.na(s)] <- ""
  s <- gsub('"', '""', s, fixed = TRUE)
  s <- gsub("[\r\n]+", " ", s)
  paste0('"', s, '"')
}

conn <- file(csvPath, open = "at")
appendRow <- function(entry, conditioning, generic, type, sample, outcome, detail) {
  writeLines(
    paste(
      escapeCsv(entry), escapeCsv(conditioning), escapeCsv(generic),
      escapeCsv(type), escapeCsv(sample), escapeCsv(outcome), escapeCsv(detail),
      sep = ","
    ),
    conn
  )
  flush(conn)
}

rawErrorPatterns <- c(
  "^'arg' should be one of", "unused argument", "subscript out of bounds",
  "could not find function", "is missing, with no default",
  "missing value where TRUE/FALSE needed", "non-numeric argument",
  "invalid subscript", "object '.*' not found", "attempt to apply non-function",
  "[$] operator is invalid for atomic vectors", "invalid 'type'", "NA/NaN/Inf",
  "unable to find an inherited method", "no applicable method",
  "incorrect number of dimensions", "non-conformable", "cannot coerce",
  "undefined columns selected", "replacement has .* rows, data has",
  "is not a function"
)
classifyError <- function(msg) {
  for (p in rawErrorPatterns) if (grepl(p, msg, perl = TRUE)) return("error-without-reason")
  "refused"
}
describeValue <- function(v) {
  if (is.null(v)) return("NULL")
  cls <- paste(class(v), collapse = "/")
  d <- tryCatch(dim(v), error = function(e) NULL)
  shape <- if (!is.null(d)) paste(d, collapse = "x") else paste0("len", length(v))
  trimws(paste(cls, shape))
}
callMethod <- function(generic, expr) {
  out <- withCallingHandlers(
    tryCatch(
      list(ok = TRUE, val = force(expr), msg = NA_character_),
      error = function(e) list(ok = FALSE, val = NULL, msg = conditionMessage(e))
    ),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  if (out$ok) {
    appendRow("dbartsSampler-method", "gaussian sampler", generic, NA, NA, "accepted", describeValue(out$val))
  } else {
    appendRow("dbartsSampler-method", "gaussian sampler", generic, NA, NA, classifyError(out$msg), out$msg)
  }
  invisible(out)
}

## ---- fixture: a plain single-forest gaussian sampler, keepTrees so trees
## and getTrees()/plotTree() have something to read ----
set.seed(20260824)
n <- 40L; nTest <- 10L
x <- matrix(runif(n * 3L), n, 3L); colnames(x) <- c("x1", "x2", "x3")
xTest <- matrix(runif(nTest * 3L), nTest, 3L); colnames(xTest) <- colnames(x)
y <- as.numeric(2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3))

ctrl <- dbartsControl(
  n.chains = 1L, n.threads = 1L, n.trees = 4L, n.burn = 3L, n.samples = 5L,
  keepTrees = TRUE, keepTrainingFits = TRUE
)
s <- dbarts(x, y, control = ctrl, verbose = FALSE)

runResult <- s$run(3L, 5L)

## ---- the 41 substantive methods, in test-host-shell-pins.R's own order ----
callMethod("run", s$run(0L, 2L))
callMethod("sampleTreesFromPrior", s$sampleTreesFromPrior())
callMethod("sampleNodeParametersFromPrior", s$sampleNodeParametersFromPrior())
callMethod("growFromRoot", s$growFromRoot(1L))
callMethod("setControl", s$setControl(ctrl))
callMethod("setModel", s$setModel(s$model))
callMethod("setData", s$setData(s$data))
callMethod("setResponse", s$setResponse(rnorm(n)))
callMethod("setOffset", s$setOffset(rep(0, n)))
callMethod("setWeights", s$setWeights(runif(n, 0.5, 1.5)))
callMethod("setCounts", s$setCounts(matrix(1L, n, 2L)))
callMethod("setCategoryOffset", s$setCategoryOffset(rep(0, n)))
callMethod("setCategoryTestOffset", s$setCategoryTestOffset(rep(0, nTest)))
callMethod("setActiveRows", s$setActiveRows(rep(TRUE, n)))
callMethod("setForestWeights", s$setForestWeights(1L, runif(n)))
callMethod("setForestBasis", s$setForestBasis(1L, ~x1))
callMethod("setSigma", s$setSigma(1.0))
callMethod("setPredictor", s$setPredictor(x))
callMethod("setCutPoints", s$setCutPoints(sort(runif(20)), column = 1L))
callMethod("setTestPredictor", { s$setTestPredictorAndOffset(xTest, NULL); s$setTestPredictor(xTest, column = 1L) })
callMethod("setTestPredictorAndOffset", s$setTestPredictorAndOffset(xTest, rep(0, nTest)))
callMethod("setTestOffset", s$setTestOffset(rep(0, nTest)))
callMethod("setCalibration", s$setCalibration(prior.scale = 1.0))
callMethod("installTrees", { s2 <- s$copy(shallow = TRUE); s$installTrees(s2) })
callMethod("getDispersion", s$getDispersion())
callMethod("getFitsWithoutOffset", s$getFitsWithoutOffset())
callMethod("copy", s$copy(shallow = TRUE))
callMethod("predict", s$predict(xTest))
callMethod("predictForests", s$predictForests(xTest))
callMethod("getLatents", s$getLatents(runResult))
callMethod("getSigmas", s$getSigmas(runResult))
callMethod("getSumsOfSquaredResiduals", s$getSumsOfSquaredResiduals(runResult))
callMethod("getForestFits", s$getForestFits(1L))
callMethod("getForestAmplitudes", s$getForestAmplitudes())
callMethod("getForestVariableCounts", s$getForestVariableCounts(1L))
callMethod("getCalibration", s$getCalibration(1L))
## setState needs an actual saved-state SEXP; storeState() first, as
## test-host-shell-pins.R's save/reload path does, then round-trip it
callMethod("storeState", s$storeState())
callMethod("setState", s$setState(s$getPointer()))
callMethod("printTrees", s$printTrees())
callMethod("getTrees", s$getTrees())
grDevices::pdf(NULL)
callMethod("plotTree", s$plotTree(1L, 1L, 1L))

close(conn)
cat("dbartsSampler method census appended to", csvPath, "\n")
