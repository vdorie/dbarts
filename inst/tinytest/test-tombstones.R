# The tombstone registry is the one list of everything kept reachable for a
# release past its removal. It has to agree with what the package actually
# exports and with what the news file tells a user, or a stub outlives its
# announcement (or, worse, is announced and not there).

registry <- dbarts:::dbartsTombstones
expect_true(length(registry) > 0L)

field <- function(name) vapply(registry, `[[`, character(1L), name)

names.t <- field("name")
kinds <- field("kind")
expires <- field("expires")

# every entry carries the four fields the registry is read by
expect_true(all(vapply(
  registry,
  function(e) {
    all(c("name", "kind", "owner", "successor", "expires") %in% names(e))
  },
  logical(1L)
)))
expect_true(all(nzchar(names.t)))
expect_true(all(
  kinds %in%
    c("function", "method", "rcMethod", "argument", "family", "behaviour")
))

# one expiry for the whole set: the release that drops them deletes the file,
# so a per-entry date would be a lie the next reader has to check
expect_identical(unique(expires), dbarts:::tombstoneExpiry)
expect_true(
  package_version(dbarts:::tombstoneExpiry) >
    package_version(as.character(packageVersion("dbarts")))
)

# --- the registry agrees with NAMESPACE ---
namespaceFile <- file.path(
  dirname(system.file("DESCRIPTION", package = "dbarts")),
  "NAMESPACE"
)
expect_true(file.exists(namespaceFile))
namespaceText <- paste(readLines(namespaceFile), collapse = "\n")

exported <- getNamespaceExports("dbarts")
for (entry in registry[kinds == "function"]) {
  expect_true(entry$name %in% exported)
  expect_true(grepl(
    paste0("\\b", entry$name, "\\b"),
    namespaceText
  ))
}

# a method tombstone is registered on its generic for its own class, and the
# registration is what makes the stub reachable at all
for (entry in registry[kinds == "method"]) {
  generic <- sub("\\.[^.]+$", "", entry$name)
  expect_true(grepl(
    paste0("S3method(", generic, ", ", entry$owner, ")"),
    namespaceText,
    fixed = TRUE
  ))
  expect_inherits(getFromNamespace(entry$name, "dbarts"), "function")
}

# a reference-class tombstone is a live method on its generator
for (entry in registry[kinds == "rcMethod"]) {
  generator <- getFromNamespace("dbartsSampler", "dbarts")
  expect_true(entry$name %in% generator$methods())
}

# an argument tombstone is reachable on its owner: either as a formal still
# carried under the old name, or through that entry point's '...'
for (entry in registry[kinds == "argument"]) {
  owner <- getFromNamespace(entry$owner, "dbarts")
  ownFormals <- names(formals(owner))
  expect_true(entry$name %in% ownFormals || "..." %in% ownFormals)
  # the successor is either a formal that still stands on this entry point
  # (a rename) or the argument of an object the setting moved onto, whose
  # own formal is named first
  target <- sub(" =.*$", "", entry$successor)
  expect_true(entry$successor %in% ownFormals || target %in% ownFormals)
}

# --- rbart_vi is live: it warns of its deprecation, then fits ---
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- NULL
rbartWarnings <- character()
rbartFit <- withCallingHandlers(
  dbarts::rbart_vi(
    y ~ x,
    data.frame(y = rnorm(30L), x = rnorm(30L)),
    group.by = rep(1:3, 10L),
    n.samples = 2L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  ),
  warning = function(w) {
    rbartWarnings <<- c(rbartWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(rbartWarnings), 1L)
expect_true(grepl("stan4bart", rbartWarnings, fixed = TRUE))
expect_true(grepl(dbarts:::tombstoneExpiry, rbartWarnings, fixed = TRUE))
printedRbart <- capture.output(print(rbartFit))
expect_true(any(grepl("rbart_vi(", printedRbart, fixed = TRUE)))
rm(onceState, rbartWarnings, rbartFit, printedRbart)

# the two thread methods are no-ops, warned once, not errors: a 0.9-x Gibbs
# loop that brackets its sweeps with them still runs
threadSampler <- dbarts::dbarts(
  matrix(rnorm(40L), 20L, 2L),
  rnorm(20L),
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 3L,
    n.samples = 2L,
    updateState = FALSE
  )
)
warnEnv <- dbarts:::onceWarnState
warnEnv[["tombstone.startThreads"]] <- NULL
warnEnv[["tombstone.stopThreads"]] <- NULL
expect_warning(threadSampler$startThreads(), pattern = "does nothing")
expect_warning(threadSampler$stopThreads(), pattern = "does nothing")
expect_null(suppressWarnings(threadSampler$startThreads()))
expect_null(suppressWarnings(threadSampler$stopThreads()))

# run's thread count, by either name or as 0.9-x's fourth positional
# argument, is ignored after a warning; anything else in its dots is refused
warnEnv[["tombstone.run.n.threads"]] <- NULL
expect_warning(
  threadSampler$run(0L, 2L, n.threads = 2L),
  pattern = "setControl"
)
expect_silent(threadSampler$run(0L, 2L, numThreads = 2L))
expect_silent(threadSampler$run(0L, 2L, NULL, 2L))
# 0.9-34's positional NA is the third argument, updateState: read as NULL,
# after its own warning
warnEnv[["tombstone.NA.updateState.dbartsSampler"]] <- NULL
expect_warning(
  threadSampler$run(0L, 2L, NA, 2L),
  pattern = "'updateState = NA' is now 'updateState = NULL'"
)
expect_error(threadSampler$run(0L, 2L, nthreads = 2L), pattern = "'nthreads'")
expect_error(threadSampler$run(0L, 2L, NULL, 2L, NULL), pattern = "<unnamed>")

# --- the consolidated argument names (dec-B98) ---

# each is accepted for one release, warned about once, and MAPPED onto the
# object it now rides, so the run is the one the old spelling asked for
consolidated <- dbarts:::onceWarnState
resetConsolidatedWarning <- function(name, caller) {
  consolidated[[paste0("tombstone.consolidated.", name, ".", caller)]] <- NULL
}

set.seed(717L)
nCons <- 40L
xCons <- matrix(rnorm(nCons * 2L), nCons, 2L)
yCons <- xCons[, 1L] + rnorm(nCons)
consControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  updateState = FALSE
)

# --- names that never reached main (dec-2026-09-24) ---

# resid.dist, breaks, max.rows, dart, levelGibbs and prior.scale only ever
# existed on the development branch; none is in 0.9-x, so there is no
# compatibility to preserve and each is simply an unknown argument now, on
# whichever door used to carry it - a plain "unused argument" refusal, not a
# silent drop. shape, which only nbinom() takes, is refused the same way
expect_error(
  dbarts::bart(xCons, yCons, resid.dist = 1, verbose = FALSE),
  pattern = "unused argument 'resid.dist'"
)
expect_error(
  dbarts::bart(xCons, yCons, shape = 1, verbose = FALSE),
  pattern = "unused argument 'shape'"
)
expect_error(
  dbarts::bart(xCons, yCons, breaks = 1, verbose = FALSE),
  pattern = "unused argument 'breaks'"
)
expect_error(
  dbarts::bart(xCons, yCons, max.rows = 1, verbose = FALSE),
  pattern = "unused argument 'max.rows'"
)
expect_error(
  dbarts::bart(xCons, yCons, dart = TRUE, verbose = FALSE),
  pattern = "unused argument 'dart'"
)
expect_error(
  dbarts::bart(xCons, yCons, levelGibbs = TRUE, verbose = FALSE),
  pattern = "unused argument 'levelGibbs'"
)
expect_error(
  dbarts::bart(xCons, yCons, prior.scale = 1, verbose = FALSE),
  pattern = "unused argument 'prior.scale'"
)
expect_error(
  dbarts::dbarts(xCons, yCons, resid.dist = 1, control = consControl),
  pattern = "unused argument 'resid.dist'"
)
expect_error(
  dbarts::dbarts(xCons, yCons, shape = 1, control = consControl),
  pattern = "unused argument 'shape'"
)
expect_error(
  dbarts::dbarts(xCons, yCons, breaks = 1, control = consControl),
  pattern = "unused argument 'breaks'"
)
expect_error(
  dbarts::dbarts(xCons, yCons, max.rows = 1, control = consControl),
  pattern = "unused argument 'max.rows'"
)
consData <- dbarts::dbartsData(xCons, yCons)
expect_error(
  dbarts::dbartsSpec(consData, resid.dist = 1, control = consControl),
  pattern = "unused argument 'resid.dist'"
)
expect_error(
  dbarts::dbartsSpec(consData, shape = 1, control = consControl),
  pattern = "unused argument 'shape'"
)
expect_error(
  dbarts::xbart(xCons, yCons, dart = TRUE, n.reps = 1L),
  pattern = "unused argument 'dart'"
)
rm(consData)

# family = "twopart" is likewise gone: refused the same way any other
# unrecognized family token is (through match.arg), rather than a named
# retirement message
expect_error(
  dbarts::bart(xCons, yCons, family = "twopart", verbose = FALSE),
  pattern = "'family' should be one of"
)
expect_error(
  dbarts::dbarts(xCons, yCons, family = "twopart", control = consControl),
  pattern = "'family' should be one of"
)

# the residual prior's three retired spellings: each warns once, each is
# applied, and each lands the same prior the family object now carries
residPriorFit <- function(...) {
  dbarts::bart(
    xCons,
    yCons,
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    keepSampler = TRUE,
    verbose = FALSE
  )
}

resetConsolidatedWarning("resid.prior", "bart")
expect_warning(
  fitResidPrior <- residPriorFit(resid.prior = dbarts::dbartsPriors$fixed(2)),
  pattern = "family = gaussian\\(sigma"
)
expect_inherits(fitResidPrior$fit$model@resid.prior, "dbartsFixedPrior")
expect_equal(fitResidPrior$fit$model@resid.prior@value, 2)
# once per session: the second call is silent and still maps
expect_silent(
  fitResidPriorAgain <- residPriorFit(
    resid.prior = dbarts::dbartsPriors$fixed(2)
  )
)
expect_equal(fitResidPriorAgain$fit$model@resid.prior@value, 2)

resetConsolidatedWarning("sigdf", "bart")
resetConsolidatedWarning("sigquant", "bart")
expect_warning(
  fitSigdf <- residPriorFit(sigdf = 5, sigquant = 0.75),
  pattern = "family = gaussian\\(sigma"
)
expect_equal(fitSigdf$fit$model@resid.prior@df, 5)
expect_equal(fitSigdf$fit$model@resid.prior@quantile, 0.75)

# the object spelling is the same fit, draw for draw. The family vocabulary
# resolves in the argument the caller writes, so these calls are spelled out
# rather than routed through the helper's '...'
fitSigmaObject <- dbarts::bart(
  xCons,
  yCons,
  family = gaussian(sigma = chisq(5, 0.75)),
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 313L,
  keepSampler = TRUE,
  verbose = FALSE
)
expect_identical(fitSigdf$yhat.train, fitSigmaObject$yhat.train)
expect_identical(fitSigdf$sigma, fitSigmaObject$sigma)

# the residual prior written both ways: agreeing spellings are one statement
# said twice and stand, disagreeing ones are refused naming both
fitAgreed <- suppressWarnings(dbarts::bart(
  xCons,
  yCons,
  family = gaussian(sigma = chisq(5, 0.75)),
  sigdf = 5,
  sigquant = 0.75,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 313L,
  keepSampler = TRUE,
  verbose = FALSE
))
expect_equal(fitAgreed$fit$model@resid.prior@df, 5)
expect_identical(fitAgreed$sigma, fitSigdf$sigma)

# both retired shorthands are named together, in one message: naming only
# the first would have a caller delete it, rerun, and hit the same refusal
# on the other
expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    yCons,
    family = gaussian(sigma = chisq(2, 0.5)),
    sigdf = 5,
    sigquant = 0.75,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    verbose = FALSE
  )),
  pattern = "'sigdf' and 'sigquant' and the family's own 'sigma' set different residual"
)
expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    yCons,
    family = gaussian(sigma = chisq(2, 0.5)),
    sigdf = 5,
    sigquant = 0.75,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    verbose = FALSE
  )),
  pattern = "Delete 'sigdf' and 'sigquant'"
)
expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    yCons,
    family = gaussian(sigma = chisq(2, 0.5)),
    resid.prior = dbarts::dbartsPriors$fixed(2),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    verbose = FALSE
  )),
  pattern = "'resid.prior' and the family's own 'sigma' set different residual"
)

# the disagreement rule compares the RESOLVED prior objects, not the
# spelling that built them: two constructor calls read back the same when
# their arguments are supplied differently (positional vs named), and a
# prebuilt object or a variable holding one agrees the same way an inline
# constructor call would. Spelled out rather than routed through
# residPriorFit's '...', as fitSigmaObject above is.
fitSameByName <- suppressWarnings(dbarts::bart(
  xCons,
  yCons,
  family = gaussian(sigma = chisq(df = 3, quant = 0.9)),
  resid.prior = dbarts::dbartsPriors$chisq(3, 0.9),
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 313L,
  keepSampler = TRUE,
  verbose = FALSE
))
expect_equal(fitSameByName$fit$model@resid.prior@df, 3)

prebuiltChisq <- dbarts::dbartsPriors$chisq(3, 0.9)
fitPrebuilt <- suppressWarnings(dbarts::bart(
  xCons,
  yCons,
  family = gaussian(sigma = dbarts::dbartsPriors$chisq(3, 0.9)),
  resid.prior = prebuiltChisq,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 313L,
  keepSampler = TRUE,
  verbose = FALSE
))
expect_equal(fitPrebuilt$fit$model@resid.prior@df, 3)

chisqVar <- prebuiltChisq
fitVariable <- suppressWarnings(dbarts::bart(
  xCons,
  yCons,
  family = gaussian(sigma = dbarts::dbartsPriors$chisq(3, 0.9)),
  resid.prior = chisqVar,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 313L,
  keepSampler = TRUE,
  verbose = FALSE
))
expect_equal(fitVariable$fit$model@resid.prior@df, 3)
rm(fitSameByName, prebuiltChisq, fitPrebuilt, chisqVar, fitVariable)

expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    yCons,
    family = gaussian(sigma = chisq(3, 0.95)),
    resid.prior = dbarts::dbartsPriors$chisq(3, 0.9),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    verbose = FALSE
  )),
  pattern = "'resid.prior' and the family's own 'sigma' set different residual"
)

# the same retirement on the other three doors: warn once, forward the value,
# and refuse a disagreement
consControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  seed = 313L,
  updateState = FALSE
)
resetConsolidatedWarning("resid.prior", "dbarts")
expect_warning(
  samplerRetired <- dbarts::dbarts(
    xCons,
    yCons,
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  ),
  pattern = "family = gaussian\\(sigma"
)
expect_equal(samplerRetired$model@resid.prior@value, 2)
expect_silent(
  samplerRetiredAgain <- dbarts::dbarts(
    xCons,
    yCons,
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  )
)
expect_equal(samplerRetiredAgain$model@resid.prior@value, 2)
expect_error(
  dbarts::dbarts(
    xCons,
    yCons,
    family = gaussian(sigma = fixed(3)),
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  ),
  pattern = "'resid.prior' and the family's own 'sigma' set different residual"
)
# a value that is not a residual prior is refused under the old spelling too
expect_error(
  dbarts::dbarts(xCons, yCons, resid.prior = NULL, control = consControl),
  pattern = "must be a residual prior specification"
)

# dbartsSpec never carried resid.prior: refused by name, naming the successor
expect_error(
  dbarts::dbartsSpec(
    dbarts::dbartsData(xCons, yCons),
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  ),
  pattern = "unused argument 'resid.prior' passed to 'dbartsSpec'.*family = gaussian\\(sigma = \\)"
)

resetConsolidatedWarning("resid.prior", "xbart")
xbartRetiredArgs <- list(
  xCons,
  yCons,
  n.reps = 1L,
  n.samples = 5L,
  n.burn = c(5L, 2L),
  n.trees = 5L,
  n.threads = 1L,
  seed = 313L
)
expect_warning(
  lossRetired <- do.call(
    dbarts::xbart,
    c(xbartRetiredArgs, list(resid.prior = dbarts::dbartsPriors$fixed(2)))
  ),
  pattern = "family = gaussian\\(sigma"
)
expect_identical(
  lossRetired,
  do.call(
    dbarts::xbart,
    c(
      xbartRetiredArgs,
      list(
        family = dbarts::dbartsFamilies$gaussian(
          sigma = dbarts::dbartsPriors$fixed(2)
        )
      )
    )
  )
)
expect_error(
  do.call(
    dbarts::xbart,
    c(
      xbartRetiredArgs,
      list(
        resid.prior = dbarts::dbartsPriors$fixed(2),
        family = dbarts::dbartsFamilies$gaussian(
          sigma = dbarts::dbartsPriors$fixed(3)
        )
      )
    )
  ),
  pattern = "'resid.prior' and the family's own 'sigma' set different residual"
)
rm(consControl, xbartRetiredArgs)

# --- xbart's own (front-door S3) ---

# a three-element n.burn was 0.9-x's per-replication burn-in; chains are
# never carried between replications now, so the element is named and gone
expect_error(
  dbarts::xbart(xCons, yCons, n.reps = 1L, n.burn = c(2L, 1L, 1L)),
  pattern = "per-replication burn-in"
)

# the prior scalars dec-B116 moved onto the prior objects and the control
for (name in dbarts:::consolidatedPriorScalars) {
  expect_false(name %in% names(formals(dbarts::bart)))
}
expect_false("proposal.probs" %in% names(formals(dbarts::dbarts)))
# 'control' is a formal of both front doors again: xbart's refusal is reversed
expect_true("control" %in% names(formals(dbarts::bart)))
expect_true("control" %in% names(formals(dbarts::xbart)))
expect_true("proposal.probs" %in% names(formals(dbarts::dbartsControl)))

# --- node.prior -> leaf.prior (dbarts), like sigma -> sigest ---

xLV <- matrix(rnorm(60L * 2L), 60L, 2L)
yLV <- xLV[, 1L] + rnorm(60L)
lvControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  seed = 919L,
  updateState = FALSE
)

warnEnv[["tombstone.node.prior.dbarts"]] <- NULL
leafWarnings <- character(0L)
fitNodePrior <- withCallingHandlers(
  dbarts::dbarts(
    xLV,
    yLV,
    control = lvControl,
    node.prior = dbarts::dbartsPriors$normal(3)
  ),
  warning = function(w) {
    leafWarnings <<- c(leafWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(leafWarnings), 1L)
fitLeafPrior <- dbarts::dbarts(
  xLV,
  yLV,
  control = lvControl,
  leaf.prior = dbarts::dbartsPriors$normal(3)
)
expect_identical(fitNodePrior$run(5L, 5L)$train, fitLeafPrior$run(5L, 5L)$train)

# node.prior = NULL is a supplied value, not an absent one (missing() is
# FALSE), so the warning still fires and the NULL still reaches parsePriors
# as leaf.prior's value; it must refuse exactly as leaf.prior = NULL does,
# not silently fall back to leaf.prior's own normal default (F1: a plain
# matchedCall$leaf.prior <- matchedCall$node.prior assignment of NULL
# deletes the element instead of setting it)
warnEnv[["tombstone.node.prior.dbarts"]] <- NULL
expect_error(
  suppressWarnings(
    dbarts::dbarts(xLV, yLV, control = lvControl, node.prior = NULL)
  ),
  pattern = "'leaf.prior' must be a leaf prior specification"
)
expect_error(
  dbarts::dbarts(xLV, yLV, control = lvControl, leaf.prior = NULL),
  pattern = "'leaf.prior' must be a leaf prior specification"
)

# supplying both spellings is an error, even when they agree
expect_error(
  dbarts::dbarts(
    xLV,
    yLV,
    control = lvControl,
    node.prior = dbarts::dbartsPriors$normal(3),
    leaf.prior = dbarts::dbartsPriors$normal(3)
  ),
  pattern = "'node.prior' and 'leaf.prior' name the same prior"
)

specDataLV <- dbarts::dbartsData(xLV, yLV)
# dbartsSpec never carried node.prior: refused by name, naming the successor
expect_error(
  dbarts::dbartsSpec(
    specDataLV,
    control = lvControl,
    node.prior = dbarts::dbartsPriors$normal(3)
  ),
  pattern = "unused argument 'node.prior' passed to 'dbartsSpec'; the leaf prior is 'leaf.prior'"
)
# ...and beside the new spelling too, rather than as a both-spellings conflict
expect_error(
  dbarts::dbartsSpec(
    specDataLV,
    control = lvControl,
    leaf.prior = dbarts::dbartsPriors$normal(3),
    node.prior = dbarts::dbartsPriors$normal(3)
  ),
  pattern = "unused argument 'node.prior' passed to 'dbartsSpec'; the leaf prior is 'leaf.prior'"
)

# bart and xbart never carried node.prior on this branch: refused by name
# rather than tombstoned (dec-B128), the message naming the successor
expect_error(
  dbarts::bart(
    xLV,
    yLV,
    node.prior = dbarts::dbartsPriors$normal(3),
    verbose = FALSE
  ),
  pattern = "unused argument 'node.prior'"
)
expect_error(
  dbarts::bart(
    xLV,
    yLV,
    node.prior = dbarts::dbartsPriors$normal(3),
    verbose = FALSE
  ),
  pattern = "leaf.prior"
)
expect_error(
  dbarts::xbart(
    xLV,
    yLV,
    node.prior = dbarts::dbartsPriors$normal(3),
    n.reps = 1L
  ),
  pattern = "unused argument 'node.prior'"
)
expect_error(
  dbarts::xbart(
    xLV,
    yLV,
    node.prior = dbarts::dbartsPriors$normal(3),
    n.reps = 1L
  ),
  pattern = "leaf.prior"
)

# --- $sampleNodeParametersFromPrior -> $sampleLeafParametersFromPrior ---
warnEnv[["tombstone.sampleNodeParametersFromPrior"]] <- NULL
freshOld <- dbarts::dbarts(xLV, yLV, control = lvControl)
freshNew <- dbarts::dbarts(xLV, yLV, control = lvControl)
sampleWarnings <- character(0L)
withCallingHandlers(
  freshOld$sampleNodeParametersFromPrior(updateState = TRUE),
  warning = function(w) {
    sampleWarnings <<- c(sampleWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(sampleWarnings), 1L)
expect_true(grepl(
  "sampleLeafParametersFromPrior",
  sampleWarnings[1L],
  fixed = TRUE
))
expect_silent(freshNew$sampleLeafParametersFromPrior(updateState = TRUE))
expect_identical(freshOld$state, freshNew$state)
# once per session: a second old-spelling call on the same sampler, with the
# key left as the first call set it, is silent and still forwards
expect_silent(freshOld$sampleNodeParametersFromPrior(updateState = FALSE))

# --- no leaks: the ordinary paths stay silent ---
countWarnings <- function(expr) {
  n <- 0L
  withCallingHandlers(
    expr,
    warning = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  n
}
expect_equal(
  countWarnings(dbarts::bart(
    xLV,
    yLV,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  0L
)
expect_equal(
  countWarnings(dbarts::bartBT(
    xLV,
    yLV,
    ntree = 5L,
    ndpost = 5L,
    nskip = 5L,
    verbose = FALSE
  )),
  0L
)
expect_equal(
  countWarnings(dbarts::xbart(
    xLV,
    yLV,
    n.reps = 1L,
    n.samples = 5L,
    n.burn = c(5L, 2L),
    n.trees = 5L,
    n.threads = 1L
  )),
  0L
)
expect_equal(countWarnings(dbarts::dbarts(xLV, yLV, control = lvControl)), 0L)
expect_equal(
  countWarnings(dbarts::dbartsSpec(specDataLV, control = lvControl)),
  0L
)
yMultiLV <- factor(sample(c("a", "b", "c"), 60L, replace = TRUE))
expect_equal(
  countWarnings(dbarts::bart(
    xLV,
    yMultiLV,
    family = "multinomial",
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  0L
)
rm(
  xLV,
  yLV,
  lvControl,
  leafWarnings,
  fitNodePrior,
  fitLeafPrior,
  specDataLV,
  freshOld,
  freshNew,
  sampleWarnings,
  countWarnings,
  yMultiLV
)

# --- the registry agrees with the news file ---
# The 1.0-0 news entry is the user-facing half of this list; the assertion
# below is the same one, read off inst/NEWS.Rd. Until the news entry lands
# (it is written with the manual, in the last slice of this arc), the file
# carries no tombstone list to check and this file stops here rather than
# pinning a section that does not exist yet.
newsFile <- file.path(
  dirname(system.file("DESCRIPTION", package = "dbarts")),
  "NEWS.Rd"
)
expect_true(file.exists(newsFile))
newsText <- paste(readLines(newsFile), collapse = "\n")
version1 <- regmatches(
  newsText,
  regexpr(
    "(?s)CHANGES IN VERSION 1\\.0-0\\}\\{.*?CHANGES IN VERSION 0",
    newsText,
    perl = TRUE
  )
)
expect_equal(length(version1), 1L)
if (!grepl("tombstone", version1, ignore.case = TRUE)) {
  exit_file(
    "NEWS tombstone list is written in the manual/news slice of this arc"
  )
}

# The check above only proves the word "tombstone" appears somewhere in a
# ~2000-line section - every registry name recurs elsewhere in there too
# (sigma, dart, control, ... are all ordinary English or reused in other
# bullets), so grepl-ing the whole section never discriminates a name
# actually dropped from the tombstone list itself. Narrow to the one
# \item that carries the list: it opens with a fixed marker sentence and
# ends where the next \item at the same indent begins.
markerStart <- "Every name and argument kept reachable"
expect_true(
  grepl(markerStart, version1, fixed = TRUE),
  info = "the tombstone list's marker sentence in inst/NEWS.Rd moved or was reworded"
)
tombstoneItem <- regmatches(
  version1,
  regexpr(
    paste0("(?s)", markerStart, ".*?(?=\\n      \\\\item )"),
    version1,
    perl = TRUE
  )
)
expect_equal(length(tombstoneItem), 1L)
expect_true(nzchar(tombstoneItem))
for (nm in unique(names.t[kinds %in% c("function", "method", "argument")])) {
  expect_true(
    grepl(nm, tombstoneItem, fixed = TRUE),
    info = paste0("'", nm, "' is missing from NEWS's tombstone-list item")
  )
}

# retired spellings resolve through forwarding wrappers to the caller's own
# values, never to a same-named global
shadowFit <- function(...) {
  dbarts::bart(
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 313L,
    keepSampler = TRUE,
    verbose = FALSE
  )
}
shadowMiddle <- function(...) shadowFit(...)
shadowRenamed <- function(power) shadowFit(xCons, yCons, power = power)
shadowData <- data.frame(y = yCons, x = xCons[, 1L])
pw <- 9
sdf <- 9
qq <- 0.1
suppressWarnings({
  fitPower <- (function() {
    pw <- 3
    shadowFit(xCons, yCons, power = pw)
  })()
  fitPowerNested <- (function() {
    pw <- 3
    shadowMiddle(xCons, yCons, power = pw)
  })()
  fitPowerRenamed <- (function() {
    power <- 3
    shadowRenamed(power)
  })()
  fitSdf <- (function() {
    sdf <- 5
    shadowFit(xCons, yCons, sigdf = sdf, sigquant = 0.75)
  })()
  fitResid <- (function() {
    qq <- 0.75
    shadowMiddle(y ~ x, shadowData, resid.prior = chisq(5, qq))
  })()
})
expect_equal(fitPower$fit$model@tree.prior@power, 3)
expect_equal(fitPowerNested$fit$model@tree.prior@power, 3)
expect_equal(fitPowerRenamed$fit$model@tree.prior@power, 3)
# a NULL power or base is the default, at bart and bartBT
suppressWarnings({
  fitDefault <- shadowFit(xCons, yCons)
  fitNullPower <- shadowFit(xCons, yCons, power = NULL)
  fitNullBase <- shadowFit(xCons, yCons, base = NULL)
  bartBTNullPower <- dbarts::bartBT(
    xCons,
    yCons,
    ntree = 5L,
    ndpost = 5L,
    nskip = 2L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE,
    power = NULL
  )
  bartBTNullBase <- dbarts::bartBT(
    xCons,
    yCons,
    ntree = 5L,
    ndpost = 5L,
    nskip = 2L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE,
    base = NULL
  )
})
expect_equal(fitNullPower$fit$model@tree.prior@power, 2)
expect_equal(fitNullPower$fit$model@tree.prior@base, 0.95)
expect_equal(fitNullBase$fit$model@tree.prior@power, 2)
expect_equal(fitNullBase$fit$model@tree.prior@base, 0.95)
expect_identical(fitNullPower$yhat.train, fitDefault$yhat.train)
expect_identical(fitNullBase$yhat.train, fitDefault$yhat.train)
expect_identical(dim(bartBTNullPower$yhat.train), c(5L, nrow(xCons)))
expect_identical(dim(bartBTNullBase$yhat.train), c(5L, nrow(xCons)))
expect_equal(fitSdf$fit$model@resid.prior@df, 5)
expect_equal(fitResid$fit$model@resid.prior@df, 5)
expect_equal(fitResid$fit$model@resid.prior@quantile, 0.75)
# split.probs: a local beats a different global of the same name through a
# wrapper, a caller-frame variable next to the vocabulary name resolves, and
# num.vars still works
probsGlobal <- c(0.5, 0.5)
suppressWarnings({
  fitProbsLocal <- (function() {
    probsGlobal <- c(0.8, 0.2)
    shadowMiddle(xCons, yCons, split.probs = probsGlobal)
  })()
  fitProbsOnlyLocal <- (function() {
    probsOnly <- c(0.8, 0.2)
    shadowMiddle(xCons, yCons, split.probs = probsOnly)
  })()
  fitProbsVocab <- shadowMiddle(xCons, yCons, split.probs = c(3, 1) / num.vars)
  fitProbsMixed <- (function() {
    wts <- c(3, 1)
    shadowMiddle(xCons, yCons, split.probs = wts / (2 * num.vars))
  })()
})
expect_equal(fitProbsLocal$fit$model@tree.prior@splitProbabilities, c(0.8, 0.2))
expect_equal(
  fitProbsOnlyLocal$fit$model@tree.prior@splitProbabilities,
  c(0.8, 0.2)
)
expect_equal(
  fitProbsVocab$fit$model@tree.prior@splitProbabilities,
  c(0.75, 0.25)
)
expect_equal(
  fitProbsMixed$fit$model@tree.prior@splitProbabilities,
  c(0.75, 0.25)
)
rm(probsGlobal, fitProbsLocal, fitProbsOnlyLocal, fitProbsVocab, fitProbsMixed)
# identical draws to the family spelling
fitFamilyResid <- (function() {
  qq <- 0.75
  shadowMiddle(
    y ~ x,
    shadowData,
    family = gaussian(sigma = dbarts::dbartsPriors$chisq(5, qq))
  )
})()
expect_equal(fitResid$sigma, fitFamilyResid$sigma)
# and warns once, using the value
resetOnce <- dbarts:::onceWarnState
resetOnce[["tombstone.consolidated.resid.prior.bart"]] <- NULL
rm(resetOnce)
expect_warning(
  shadowFit(y ~ x, shadowData, resid.prior = chisq(3, 0.9)),
  "resid.prior"
)
rm(
  shadowFit,
  shadowMiddle,
  shadowRenamed,
  shadowData,
  pw,
  sdf,
  qq,
  fitPower,
  fitPowerNested,
  fitPowerRenamed,
  fitSdf,
  fitResid,
  fitFamilyResid
)
