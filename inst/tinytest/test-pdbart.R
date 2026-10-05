source(system.file("common", "pdData.R", package = "dbarts"), local = TRUE)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

onceState <- dbarts:::onceWarnState
# the once-per-session keys these tests clear, restored at the end of the file
savedKeys <- ls(onceState, all.names = TRUE)
savedState <- mget(savedKeys, envir = onceState)
resetPdbartKeys <- function() {
  for (key in grep("^tombstone\\.pd2?bart", ls(onceState), value = TRUE)) {
    onceState[[key]] <- NULL
  }
}

x <- testData$x
y <- testData$y
levs <- list(seq(-1, 1, 0.2), seq(-1, 1, 0.2))
small <- function(f, ...) {
  f(
    ...,
    n.trees = 5L,
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 2L,
    n.threads = 1L,
    verbose = FALSE
  )
}

# a data call, a fit with its trees, and the refit of a fit kept without them
# are one model at one seed
pdb1 <- suppressMessages(small(
  dbarts::pdbart,
  x,
  y,
  xind = c(1, 2),
  levs = levs,
  pl = FALSE,
  seed = 3L
))
expect_equal(dim(pdb1$fd[[1L]]), c(20L, 11L))
expect_identical(pdb1$n.chains, 2L)
bartFit <- small(dbarts::bart, x, y, seed = 3L, keepTrees = TRUE)
pdb2 <- dbarts::pdbart(bartFit, xind = c(1, 2), levs = levs, pl = FALSE)
expect_identical(pdb1$fd, pdb2$fd)
expect_identical(pdb2$bartcall, bartFit$call)
bartFit <- small(dbarts::bart, x, y, seed = 3L)
warnings.refit <- captureWarnings(
  pdb3 <- dbarts::pdbart(bartFit, xind = c(1, 2), levs = levs, pl = FALSE)
)
expect_equal(length(warnings.refit), 1L)
expect_inherits(warnings.refit[[1L]], "dbartsFallbackWarning")
expect_identical(pdb1$fd, pdb3$fd)
expect_identical(pdb1$bartcall[[1L]], quote(dbarts::bart))

# and under a grown initial forest
pdbGrown <- small(
  dbarts::pdbart,
  x,
  y,
  xind = 1L,
  pl = FALSE,
  seed = 3L,
  n.grow.sweeps = 2L
)
grownFit <- small(dbarts::bart, x, y, seed = 3L, n.grow.sweeps = 2L)
expect_identical(
  pdbGrown$fd,
  suppressWarnings(dbarts::pdbart(grownFit, xind = 1L, pl = FALSE))$fd
)

# each value is the per-draw mean of predictions over the fit's rows
atValue <- x
atValue[, 1L] <- levs[[1L]][4L]
fitWithTrees <- small(dbarts::bart, x, y, seed = 3L, keepTrees = TRUE)
expect_equal(pdb1$fd[[1L]][, 4L], rowMeans(predict(fitWithTrees, atValue)))

# a sampler passed in, with saved trees, against its own $predict
control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE
)
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
pdbSampler <- dbarts::pdbart(sampler, xind = 1L, levs = levs[1L], pl = FALSE)
expect_equal(pdbSampler$fd[[1L]][, 4L], colMeans(sampler$predict(atValue)))
# without saved trees it runs, with a warning, and averages its test fits
control@keepTrees <- FALSE
sampler <- dbarts::dbarts(x, y, control = control)
invisible(sampler$run())
warnings.sampler <- captureWarnings(
  pdbRun <- dbarts::pdbart(sampler, xind = 1L, levs = levs[1L], pl = FALSE)
)
expect_equal(length(warnings.sampler), 1L)
expect_inherits(warnings.sampler[[1L]], "dbartsFallbackWarning")
expect_equal(dim(pdbRun$fd[[1L]]), c(10L, 11L))
expect_false(is.null(pdbRun$yhat.train))
rm(pdbGrown, grownFit, atValue, fitWithTrees, control, sampler, pdbSampler)
# and its rows' offsets enter the averages it runs: the same sampler, seeded,
# run by hand over the same rows
offsetSampler <- function() {
  s <- dbarts::dbarts(
    x,
    y,
    offset = x[, 3L],
    control = dbarts::dbartsControl(
      n.trees = 5L,
      n.samples = 10L,
      n.burn = 5L,
      n.chains = 1L,
      n.threads = 1L,
      seed = 7L
    )
  )
  invisible(s$run())
  s
}
pdbRunOffset <- suppressWarnings(dbarts::pdbart(
  offsetSampler(),
  xind = 1L,
  levs = list(0.5),
  pl = FALSE
))
handSampler <- offsetSampler()
atHalf <- x
atHalf[, 1L] <- 0.5
handSampler$setTestPredictor(atHalf)
handSamples <- handSampler$run(0L, 10L)
expect_equal(
  pdbRunOffset$fd[[1L]][, 1L],
  colMeans(handSamples$test) + mean(x[, 3L])
)
rm(offsetSampler, pdbRunOffset, handSampler, atHalf, handSamples)

# a bartBT fit kept without trees is refit through bartBT, not forwarded to it
# through bart, which would warn a second time
onceState[["tombstone.bartShim"]] <- NULL
legacyFit <- function(...) {
  dbarts::bartBT(
    x,
    y,
    ntree = 5L,
    ndpost = 10L,
    nskip = 5L,
    seed = 3L,
    verbose = FALSE,
    ...
  )
}
warnings.bartBT <- captureWarnings(
  pdbBT <- dbarts::pdbart(legacyFit(), xind = 1L, levs = levs[1L], pl = FALSE)
)
expect_equal(length(warnings.bartBT), 1L)
expect_inherits(warnings.bartBT[[1L]], "dbartsFallbackWarning")
expect_identical(
  pdbBT$fd,
  dbarts::pdbart(
    legacyFit(keeptrees = TRUE),
    xind = 1L,
    levs = levs[1L],
    pl = FALSE
  )$fd
)
rm(legacyFit, warnings.bartBT, pdbBT)
rm(pdbRun, warnings.sampler, warnings.refit, pdb2, pdb3, bartFit)

# a fit passed in takes no fitting argument, and named x.train it is used
# with a warning to pass it first
expect_error(
  dbarts::pdbart(pdb1$fit, n.trees = 3L, pl = FALSE),
  "'n.trees' has no effect on a sampler"
)
bartFit <- small(dbarts::bart, x, y, seed = 3L, keepTrees = TRUE)
warnings.named <- captureWarnings(
  pdbNamed <- dbarts::pdbart(x.train = bartFit, xind = 1L, pl = FALSE)
)
expect_equal(length(warnings.named), 1L)
expect_true(grepl("first, unnamed", conditionMessage(warnings.named[[1L]])))
rm(warnings.named, pdbNamed)

# --- the result's shapes ---
# combineChains = FALSE gives each fd a leading chain margin, as the fit's own
# components have, and merging it back gives the default
pdbSplit <- small(
  dbarts::pdbart,
  x,
  y,
  xind = c(1, 2),
  levs = levs,
  pl = FALSE,
  seed = 3L,
  combineChains = FALSE
)
expect_equal(dim(pdbSplit$fd[[2L]]), c(2L, 10L, 11L))
expect_equal(dim(pdbSplit$yhat.train), c(2L, 10L, 100L))
expect_identical(dbarts:::combineChains(pdbSplit$fd[[2L]]), pdb1$fd[[2L]])
pd2Merged <- small(dbarts::pd2bart, x, y, xind = 2:3, pl = FALSE, seed = 3L)
pd2Split <- small(
  dbarts::pd2bart,
  x,
  y,
  xind = 2:3,
  pl = FALSE,
  seed = 3L,
  combineChains = FALSE
)
expect_equal(dim(pd2Merged$fd), c(20L, 121L))
expect_equal(dim(pd2Split$fd), c(2L, 10L, 121L))
expect_identical(dbarts:::combineChains(pd2Split$fd), pd2Merged$fd)

# keepSampler = FALSE drops the sampler, under either spelling
pdbNoSampler <- small(
  dbarts::pdbart,
  x,
  y,
  xind = 1L,
  pl = FALSE,
  seed = 3L,
  keepSampler = FALSE
)
expect_false("fit" %in% names(pdbNoSampler))
expect_false(
  "fit" %in%
    names(dbarts::pdbart(bartFit, xind = 1L, pl = FALSE, keepSampler = FALSE))
)
resetPdbartKeys()
warnings.keep <- captureWarnings(
  pdbKeep <- dbarts::pdbart(bartFit, xind = 1L, pl = FALSE, keepsampler = FALSE)
)
expect_equal(length(warnings.keep), 1L)
expect_true(grepl("'keepsampler'", conditionMessage(warnings.keep[[1L]])))
expect_false(grepl("defaults", conditionMessage(warnings.keep[[1L]])))
expect_false("fit" %in% names(pdbKeep))
rm(warnings.keep, pdbKeep)

# --- the plot methods, on both shapes ---
# one device page per entry of xind, each carrying that predictor's label
psFile <- tempfile(fileext = ".ps")
postscript(psFile, onefile = TRUE)
expect_silent(plot(pdb1))
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("^%%Page:", psLines)), 2L)
labelHits <- vapply(
  pdb1$xlbs,
  function(lab) length(grep(paste0("(", lab, ")"), psLines, fixed = TRUE)),
  integer(1L)
)
expect_equivalent(labelHits, c(1L, 1L))
# and xind selects WHICH predictor, not just how many: the second alone
postscript(psFile, onefile = TRUE)
expect_silent(plot(pdb1, xind = 2L))
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("^%%Page:", psLines)), 1L)
expect_equal(
  length(grep(paste0("(", pdb1$xlbs[2L], ")"), psLines, fixed = TRUE)),
  1L
)
# a caller's type, xlab and ylab reach the lines and the axes
postscript(psFile, onefile = TRUE)
expect_silent(plot(pdb1, xind = 1L, type = "l", xlab = "a", ylab = "b"))
dev.off()
psLines <- readLines(psFile, warn = FALSE)
expect_equal(length(grep("(a)", psLines, fixed = TRUE)), 1L)
expect_equal(length(grep("(b)", psLines, fixed = TRUE)), 1L)
expect_equal(
  length(grep(paste0("(", pdb1$xlbs[1L], ")"), psLines, fixed = TRUE)),
  0L
)
# and the type drawn is the caller's: a frame with lines only differs from
# one with lines and points
drawing <- function(...) {
  postscript(psFile, onefile = TRUE)
  plot(pdb1, xind = 1L, ...)
  dev.off()
  grep("^%%", readLines(psFile, warn = FALSE), value = TRUE, invert = TRUE)
}
expect_false(identical(drawing(type = "l"), drawing()))
expect_identical(drawing(type = "b"), drawing())
unlink(psFile)
rm(psFile, psLines, labelHits, drawing)

pdf(NULL)
expect_silent(plot(pdbSplit))
expect_silent(plot(pd2Merged))
expect_silent(plot(pd2Split, justmedian = FALSE))
dev.off()
# plot.pd2bart hands a type in its dots to image like any other argument
pdf(NULL)
imageError <- tryCatch(
  image(1:2, 1:2, matrix(1:4, 2L), type = "l"),
  error = conditionMessage
)
expect_error(plot(pd2Merged, type = "l"), imageError, fixed = TRUE)
dev.off()
rm(pdbSplit, pd2Merged, pd2Split, pdbNoSampler, imageError)

# --- argument handling ---
# pdbart's own settings are refused under bart's names and BayesTree's
for (name in c("samplerOnly", "sampleronly", "test", "x.test", "offset.test")) {
  args <- list(x, y, pl = FALSE, TRUE)
  names(args)[4L] <- name
  expect_error(do.call(dbarts::pdbart, args), "set internally")
}
expect_error(
  dbarts::pdbart(x, y, keepTrees = FALSE, pl = FALSE),
  "'keepTrees' can only be TRUE"
)
expect_error(
  dbarts::pd2bart(x, y, keeptrees = FALSE, pl = FALSE),
  "'keeptrees' can only be TRUE"
)
# a misspelled argument is refused where it used to run the default fit
expect_error(dbarts::pdbart(x, y, ntrees = 5L, pl = FALSE), "'ntrees'")
# both spellings of one setting
expect_error(
  dbarts::pdbart(x, y, ntree = 5L, n.trees = 5L, pl = FALSE),
  "'ntree' and 'n.trees'"
)
expect_error(
  dbarts::pdbart(x, y, power = 3, tree.prior = cgm(), pl = FALSE),
  "'power' and 'tree.prior'"
)
# n.trees and family are honoured
pdbProbit <- dbarts::pdbart(
  x,
  as.numeric(y > 0),
  xind = 1L,
  pl = FALSE,
  n.trees = 7L,
  n.samples = 2L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  family = "logistic",
  verbose = FALSE
)
expect_identical(pdbProbit$fit$control@n.trees, 7L)
expect_identical(pdbProbit$fit$model@family, "logistic")
rm(pdbProbit, name, args)

# The BayesTree spellings: the table is every name bartBT takes and bart does
# not, less the two pdbart sets itself
expect_true(setequal(
  names(dbarts:::pdbartBayesTreeNames),
  setdiff(
    setdiff(names(formals(dbarts::bartBT)), names(formals(dbarts::bart))),
    c("x.test", "sampleronly")
  )
))
expect_true(all(
  dbarts:::pdbartBayesTreeNames %in%
    c(names(formals(dbarts::bart)), "sigdf", "sigquant")
))

# each is translated to the setting it names, with one warning per name and
# none from bart
resetPdbartKeys()
consolidatedKey <- "tombstone.consolidated.sigdf.bart"
consolidatedBefore <- onceState[[consolidatedKey]]
legacyArgs <- list(
  ntree = 3L,
  ndpost = 4L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  printevery = 7L,
  keepevery = 2L,
  numcut = 13L,
  usequants = TRUE,
  keeptrainfits = FALSE,
  printcutoffs = 0L,
  combinechains = TRUE,
  keeptrees = TRUE,
  keepcall = TRUE,
  keepsampler = TRUE,
  sigdf = 5,
  sigquant = 0.5,
  power = 3,
  base = 0.5,
  splitprobs = c(0.7, 0.1, 0.2),
  proposalprobs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
)
warnings.legacy <- captureWarnings(
  pdbLegacy <- do.call(
    dbarts::pdbart,
    c(
      list(x.train = x, y.train = y, xind = 1L, pl = FALSE, verbose = FALSE),
      legacyArgs
    )
  )
)
expect_equal(length(warnings.legacy), length(legacyArgs) + 2L)
expect_true(all(vapply(
  warnings.legacy,
  inherits,
  FALSE,
  "dbartsDeprecatedWarning"
)))
expect_identical(onceState[[consolidatedKey]], consolidatedBefore)
legacyControl <- pdbLegacy$fit$control
expect_identical(legacyControl@n.trees, 3L)
expect_identical(legacyControl@n.chains, 1L)
expect_identical(legacyControl@n.threads, 1L)
expect_identical(legacyControl@n.thin, 2L)
expect_identical(legacyControl@n.samples, 2L)
expect_identical(legacyControl@printEvery, 3L)
expect_identical(legacyControl@n.cuts, 13L)
expect_true(legacyControl@useQuantiles)
expect_false(legacyControl@keepTrainingFits)
expect_null(pdbLegacy$yhat.train)
expect_false(is.null(pdbLegacy$fit))
expect_equal(
  legacyControl@proposal.probs[c("birth_death", "swap", "change")],
  c(birth_death = 0.5, swap = 0.1, change = 0.4)
)
expect_equal(pdbLegacy$fit$model@tree.prior@power, 3)
expect_equal(pdbLegacy$fit$model@tree.prior@base, 0.5)
expect_equal(
  pdbLegacy$fit$model@tree.prior@splitProbabilities,
  c(0.7, 0.1, 0.2)
)
expect_equal(pdbLegacy$fit$model@resid.prior@df, 5)
expect_equal(pdbLegacy$fit$model@resid.prior@quantile, 0.5)
expect_identical(pdbLegacy$bartcall$n.trees, 3L)
expect_true(all(
  c("formula", "data", "tree.prior") %in% names(pdbLegacy$bartcall)
))
expect_false(any(
  setdiff(names(legacyArgs), c("sigdf", "sigquant")) %in%
    names(pdbLegacy$bartcall)
))
# a second call in the session is silent
expect_equal(
  length(captureWarnings(do.call(
    dbarts::pdbart,
    c(list(x, y, xind = 1L, pl = FALSE, verbose = FALSE), legacyArgs)
  ))),
  0L
)
# keepcall = FALSE and keepsampler = FALSE
resetPdbartKeys()
pdbBare <- suppressWarnings(dbarts::pd2bart(
  x,
  y,
  pl = FALSE,
  ntree = 3L,
  ndpost = 2L,
  nskip = 0L,
  nchain = 1L,
  nthread = 1L,
  keepcall = FALSE,
  keepsampler = FALSE,
  verbose = FALSE
))
expect_identical(pdbBare$bartcall, call("NULL"))
expect_false("fit" %in% names(pdbBare))
# proposalprobs is set on a copy of a given control, its other settings
# standing
resetPdbartKeys()
pdbControl <- suppressWarnings(dbarts::pdbart(
  x,
  y,
  xind = 1L,
  pl = FALSE,
  ntree = 3L,
  ndpost = 2L,
  nskip = 0L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  control = dbarts::dbartsControl(n.cuts = 9L),
  proposalprobs = c(birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5)
))
expect_identical(pdbControl$fit$control@n.cuts, 9L)
expect_equal(pdbControl$fit$control@proposal.probs[["birth_death"]], 0.6)
# binaryOffset is the offset
resetPdbartKeys()
pdbOffset <- suppressWarnings(dbarts::pdbart(
  x,
  as.numeric(y > 0),
  xind = 1L,
  pl = FALSE,
  ntree = 3L,
  ndpost = 2L,
  nskip = 0L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  binaryOffset = 2
))
atValue <- x
atValue[, 1L] <- pdbOffset$levs[[1L]][2L]
expect_equal(
  pdbOffset$fd[[1L]][, 2L],
  colMeans(pdbOffset$fit$predict(atValue)) + 2
)
rm(pdbLegacy, legacyControl, legacyArgs, warnings.legacy, pdbBare)
rm(pdbControl, pdbOffset, atValue, consolidatedKey, consolidatedBefore)

# the translation warning reaches a call made from package code
resetPdbartKeys()
packageCaller <- function(x, y) {
  dbarts::pdbart(
    x,
    y,
    xind = 1L,
    pl = FALSE,
    ntree = 3L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}
environment(packageCaller) <- asNamespace("stats")
warnings.package <- captureWarnings(packageCaller(x, y))
expect_equal(length(warnings.package), 1L)
expect_true(grepl("'ntree'", conditionMessage(warnings.package[[1L]])))
rm(packageCaller, warnings.package)

# an argument the caller wrote unevaluated reaches bart as written: a local
# 'chi' stub, as treatSens keeps, does not shadow the leaf prior's
chiCaller <- function(x, y) {
  chi <- function(df, scale) stop("the stub was called")
  dbarts::pdbart(
    x,
    y,
    xind = 1L,
    pl = FALSE,
    k = chi(1.25, Inf),
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}
pdbChi <- chiCaller(x, as.numeric(y > 0))
leafPrior <- pdbChi$fit$getLeafPrior()$leaf.prior
expect_inherits(leafPrior@k, "dbartsChiHyperprior")
expect_equal(leafPrior@k@degreesOfFreedom, 1.25)
expect_equal(leafPrior@k@scale, Inf)
rm(chiCaller, pdbChi, leafPrior)

# names pdbart does not translate still warn through bart: only the notices
# pdbart replaces are held back
for (name in c("resid.prior", "split.probs")) {
  onceState[[paste0("tombstone.consolidated.", name, ".bart")]] <- NULL
}
warnings.retired <- captureWarnings(suppressMessages(dbarts::pdbart(
  x,
  y,
  xind = 1L,
  pl = FALSE,
  resid.prior = chisq(3, 0.9),
  split.probs = c(0.5, 0.25, 0.25),
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)))
expect_equal(length(warnings.retired), 2L)
retiredMessages <- vapply(warnings.retired, conditionMessage, "")
expect_equal(sum(grepl("'resid.prior'", retiredMessages)), 1L)
expect_equal(sum(grepl("'split.probs'", retiredMessages)), 1L)
rm(name, warnings.retired, retiredMessages)

# --- the defaults message ---
# once per session, for a data call naming no BayesTree argument, and not for
# package code; bart's own message is held back without using up its showing
messageKeys <- c(
  dbarts:::pdbartDefaultsKey,
  dbarts:::frontDoorDefaultsKey
)
for (key in messageKeys) {
  onceState[[key]] <- NULL
}
countMessages <- function(expr) {
  observed <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      observed <<- c(observed, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  observed
}
# a BayesTree-spelled call has its translation warning instead
expect_equal(
  length(countMessages(suppressWarnings(dbarts::pdbart(
    x,
    y,
    xind = 1L,
    pl = FALSE,
    ntree = 3L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )))),
  0L
)
expect_false(isTRUE(onceState[[dbarts:::pdbartDefaultsKey]]))
quietCaller <- function(x, y) {
  dbarts::pdbart(
    x,
    y,
    xind = 1L,
    pl = FALSE,
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}
environment(quietCaller) <- asNamespace("stats")
expect_equal(length(countMessages(quietCaller(x, y))), 0L)
messages.first <- countMessages(small(
  dbarts::pdbart,
  x,
  y,
  xind = 1L,
  pl = FALSE
))
expect_equal(length(messages.first), 1L)
expect_true(grepl("'pdbart' and 'pd2bart'", messages.first))
expect_equal(
  length(countMessages(small(dbarts::pd2bart, x, y, pl = FALSE))),
  0L
)
messages.bart <- countMessages(small(dbarts::bart, x, y))
expect_equal(length(messages.bart), 1L)
expect_true(grepl("'bart2'", messages.bart))
rm(key, messageKeys, countMessages, quietCaller, messages.first)
rm(messages.bart)

# --- families ---
# every family refusal comes before anything is fitted: a fit would refuse
# n.trees = -1 first
noFit <- function(...) {
  dbarts::pdbart(..., pl = FALSE, n.trees = -1L, verbose = FALSE)
}
set.seed(5)
n <- nrow(x)
counts <- matrix(rpois(3L * n, 2), n, 3L)
expect_error(
  noFit(x, factor(sample(letters[1:3], n, TRUE))),
  "multinomial fit"
)
expect_error(
  noFit(x, factor(sample(letters[1:3], n, TRUE), ordered = TRUE)),
  "ordinal fit"
)
expect_error(noFit(x, counts), "multinomial fit")
expect_error(
  dbarts::pd2bart(
    x,
    y,
    family = "multinomial",
    pl = FALSE,
    n.trees = -1L,
    verbose = FALSE
  ),
  "'pd2bart' does not serve a multinomial fit"
)

# and a fit or sampler of each refused family passed in
tiny <- function(response, ...) {
  dbarts::bart(
    x,
    response,
    ...,
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
refusedFits <- list(
  multinomial = tiny(factor(sample(letters[1:3], n, TRUE))),
  ordinal = tiny(factor(sample(letters[1:3], n, TRUE), ordered = TRUE))
)
for (family in names(refusedFits)) {
  refusedFit <- refusedFits[[family]]
  expect_error(
    dbarts::pdbart(refusedFit, pl = FALSE),
    "does not serve",
    info = family
  )
  if (!is.null(refusedFit$fit)) {
    expect_error(
      dbarts::pd2bart(refusedFit$fit, pl = FALSE),
      "does not serve",
      info = family
    )
  }
}
rm(noFit, n, counts, tiny, refusedFits, family, refusedFit)

# a training column with a missing value still gives a grid
xMissing <- x
xMissing[c(3L, 17L), 1L] <- NA
pdbMissing <- small(dbarts::pdbart, xMissing, y, xind = 1L, pl = FALSE)
expect_equal(
  pdbMissing$levs[[1L]],
  quantile(xMissing[, 1L], c(0.05, seq(0.1, 0.9, 0.1), 0.95), na.rm = TRUE),
  check.attributes = FALSE
)
expect_false(anyNA(pdbMissing$fd[[1L]]))
rm(xMissing, pdbMissing)

# --- offsets and the averaged rows ---
# each row's offset is added before averaging, as predict adds the fit's own
offsetFit <- small(
  dbarts::bart,
  x,
  y,
  offset = x[, 3L],
  seed = 3L,
  keepTrees = TRUE
)
pdbOff <- dbarts::pdbart(offsetFit, xind = 1L, levs = levs[1L], pl = FALSE)
atValue <- x
atValue[, 1L] <- levs[[1L]][2L]
expect_equal(pdbOff$fd[[1L]][, 2L], rowMeans(predict(offsetFit, atValue)))
expect_equal(
  pdbOff$fd[[1L]][, 2L],
  as.vector(colMeans(offsetFit$fit$predict(atValue))) + mean(x[, 3L])
)
# and on a probit fit, on the link scale
yBinary <- as.numeric(y > 0)
probitFit <- small(
  dbarts::bart,
  x,
  yBinary,
  offset = rep(2, nrow(x)),
  seed = 3L,
  keepTrees = TRUE
)
pdbProbit <- dbarts::pdbart(probitFit, xind = 1L, levs = levs[1L], pl = FALSE)
expect_equal(
  pdbProbit$fd[[1L]][, 2L],
  rowMeans(predict(probitFit, atValue, type = "bart"))
)
expect_equal(
  pdbProbit$fd[[1L]][, 2L],
  as.vector(colMeans(probitFit$fit$predict(atValue))) + 2
)
# a probit fit's 0/1 weights mask rows out, and those are not averaged
active <- rep_len(c(1, 0, 1), nrow(x))
maskedFit <- small(
  dbarts::bart,
  x,
  yBinary,
  weights = active,
  seed = 3L,
  keepTrees = TRUE
)
pdbMasked <- dbarts::pdbart(maskedFit, xind = 1L, levs = levs[1L], pl = FALSE)
expect_equal(
  pdbMasked$fd[[1L]][, 2L],
  rowMeans(predict(maskedFit, atValue[active == 1, ], type = "bart"))
)
# and a gaussian fit's 0 weights likewise
weightedFit <- withCallingHandlers(
  small(dbarts::bart, x, y, weights = active * 2, seed = 3L, keepTrees = TRUE),
  warning = function(w) {
    if (grepl("'weights' of 0", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
)
pd2Weighted <- dbarts::pd2bart(
  weightedFit,
  xind = 1:2,
  levs = list(c(-0.5, 0.5), c(0, 0.5)),
  pl = FALSE
)
atPoint <- x[active == 1, ]
atPoint[, 1L] <- 0.5
atPoint[, 2L] <- 0
expect_equal(pd2Weighted$fd[, 2L], rowMeans(predict(weightedFit, atPoint)))
rm(offsetFit, pdbOff, atValue, yBinary, probitFit, pdbProbit, active)
rm(maskedFit, pdbMasked, weightedFit, pd2Weighted, atPoint)

# --- pd2bart ---
pd2b1 <- small(
  dbarts::pd2bart,
  x,
  y,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  seed = 3L
)
bartFit <- small(dbarts::bart, x, y, seed = 3L)
pd2b2 <- suppressWarnings(dbarts::pd2bart(
  bartFit,
  xind = c(2, 3),
  pl = FALSE,
  levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
))
expect_identical(pd2b1$fd, pd2b2$fd)
pdf(NULL)
expect_silent(plot(pd2b1))
dev.off()
rm(pd2b1, pd2b2, bartFit)

# A factor predictor is evaluated at every level by default, levels are given
# and reported by name, and each value equals a direct prediction at the level.
set.seed(7)
factorFrame <- data.frame(
  y = rnorm(150L),
  x1 = runif(150L),
  f = factor(sample(LETTERS[1:15], 150L, TRUE))
)
factorFrame$y <- factorFrame$y + as.integer(factorFrame$f)
factorFit <- bart(
  y ~ .,
  factorFrame,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
pd <- pdbart(factorFit, xind = "f", pl = FALSE)
expect_identical(pd$levs[[1L]], LETTERS[1:15])
atLevel <- factorFrame
atLevel$f <- factor("C", levels = levels(factorFrame$f))
expect_equal(pd$fd[[1L]][, 3L], rowMeans(predict(factorFit, atLevel)))
pd <- pdbart(factorFit, xind = "f", levs = list(c("D", "A")), pl = FALSE)
expect_identical(pd$levs[[1L]], c("D", "A"))
atLevel$f <- factor("A", levels = levels(factorFrame$f))
expect_equal(pd$fd[[1L]][, 2L], rowMeans(predict(factorFit, atLevel)))
expect_error(
  pdbart(factorFit, xind = "f", levs = list(c("A", "ZZ")), pl = FALSE),
  "'ZZ'"
)
expect_error(
  pdbart(factorFit, xind = "f", levs = list(1:2), pl = FALSE),
  "must name its levels"
)
pd2 <- pd2bart(factorFit, xind = c("x1", "f"), pl = FALSE)
expect_identical(pd2$levs[[2L]], LETTERS[1:15])
expect_equal(ncol(pd2$fd), length(pd2$levs[[1L]]) * 15L)
grDevices::pdf(NULL)
plot(pd)
plot(pd2)
grDevices::dev.off()
rm(factorFrame, factorFit, pd, atLevel, pd2)

# With exactly two predictors and no offset each grid point is a whole row:
# the result keeps one column per grid point, equal to a direct prediction
# there, for named predictors, a fit object, and a sampler with or without
# saved trees.
set.seed(8)
xTwo <- matrix(runif(100L), 50L, 2L, dimnames = list(NULL, c("x1", "x2")))
yTwo <- 2 * xTwo[, 1L] + rnorm(50L, sd = 0.1)
twoControl <- dbarts::dbartsControl(
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L
)
for (keepTrees in c(TRUE, FALSE)) {
  twoControl@keepTrees <- keepTrees
  twoSampler <- dbarts::dbarts(xTwo, yTwo, control = twoControl)
  invisible(twoSampler$run())
  pd2 <- suppressWarnings(pd2bart(twoSampler, xind = c(1L, 2L), pl = FALSE))
  expect_equal(dim(pd2$fd), c(10L, 121L))
}
grDevices::pdf(NULL)
plot(pd2)
grDevices::dev.off()
twoFit <- bart(
  xTwo,
  yTwo,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
pd2 <- pd2bart(twoFit, xind = c(2L, 1L), pl = FALSE)
grid <- as.matrix(expand.grid(pd2$levs[[1L]], pd2$levs[[2L]]))[, c(2L, 1L)]
colnames(grid) <- c("x1", "x2")
expect_equal(pd2$fd, predict(twoFit, grid), check.attributes = FALSE)
# a varying offset makes the rows differ, so the average is taken over them
twoOffset <- bart(
  xTwo,
  yTwo,
  offset = xTwo[, 1L],
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
pd2 <- pd2bart(
  twoOffset,
  xind = 1:2,
  levs = list(c(0.2, 0.4), c(0.1, 0.9)),
  pl = FALSE
)
atPoint <- xTwo
atPoint[, 1L] <- 0.4
atPoint[, 2L] <- 0.9
expect_equal(pd2$fd[, 4L], rowMeans(predict(twoOffset, atPoint)))
rm(xTwo, yTwo, twoControl, keepTrees, twoSampler, pd2, twoFit, grid)
rm(twoOffset, atPoint)

rm(pdb1, x, y, levs, small, testData)
for (key in setdiff(ls(onceState, all.names = TRUE), savedKeys)) {
  onceState[[key]] <- NULL
}
for (key in savedKeys) {
  onceState[[key]] <- savedState[[key]]
}
rm(key, onceState, savedKeys, savedState, resetPdbartKeys)
