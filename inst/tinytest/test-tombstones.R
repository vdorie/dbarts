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

# --- the stubs themselves ---
expect_error(dbarts::rbart_vi(), pattern = "stan4bart")
expect_error(dbarts::rbart_vi(), pattern = "results move")
fakeRbart <- structure(list(), class = "rbart")
expect_error(predict(fakeRbart), pattern = "stan4bart")
expect_error(fitted(fakeRbart), pattern = "stan4bart")
expect_error(residuals(fakeRbart), pattern = "stan4bart")
expect_error(dbarts::extract(fakeRbart), pattern = "stan4bart")

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

resetConsolidatedWarning("resid.dist", "dbarts")
expect_warning(
  samplerResidDist <- dbarts::dbarts(
    xCons,
    yCons,
    control = consControl,
    resid.dist = student(df = 5)
  ),
  pattern = "family = student"
)
expect_equal(attr(samplerResidDist$model, "resid.df"), 5)
# once per session: the second call is silent and still maps
expect_silent(
  samplerResidDistAgain <- dbarts::dbarts(
    xCons,
    yCons,
    control = consControl,
    resid.dist = student(df = 5)
  )
)
expect_equal(attr(samplerResidDistAgain$model, "resid.df"), 5)

resetConsolidatedWarning("dart", "bart")
expect_warning(
  fitDart <- dbarts::bart(
    xCons,
    yCons,
    dart = TRUE,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    keepSampler = TRUE,
    verbose = FALSE
  ),
  pattern = "tree.prior = dart"
)
expect_inherits(fitDart$fit$model@tree.prior, "dbartsDartPrior")

resetConsolidatedWarning("levelGibbs", "bart")
expect_warning(
  fitLevelGibbs <- dbarts::bart(
    xCons,
    yCons,
    levelGibbs = TRUE,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    keepSampler = TRUE,
    verbose = FALSE
  ),
  pattern = "tree.prior = cgm"
)
expect_true(fitLevelGibbs$fit$control@levelGibbs)

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

resetConsolidatedWarning("resid.prior", "dbartsSpec")
expect_warning(
  specRetired <- dbarts::dbartsSpec(
    dbarts::dbartsData(xCons, yCons),
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  ),
  pattern = "family = gaussian\\(sigma"
)
expect_equal(specRetired$model@resid.prior@value, 2)
expect_error(
  dbarts::dbartsSpec(
    dbarts::dbartsData(xCons, yCons),
    family = gaussian(sigma = fixed(3)),
    resid.prior = dbarts::dbartsPriors$fixed(2),
    control = consControl
  ),
  pattern = "'resid.prior' and the family's own 'sigma' set different residual"
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

# an old spelling that names a family the call cannot fit is refused rather
# than resolved one way in silence
resetConsolidatedWarning("resid.dist", "bart")
expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    as.numeric(yCons > 0),
    family = "probit",
    resid.dist = student(df = 5),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  pattern = "student residuals require a continuous gaussian response"
)

# --- xbart's own (front-door S3) ---

# a three-element n.burn was 0.9-x's per-replication burn-in; chains are
# never carried between replications now, so the element is named and gone
expect_error(
  dbarts::xbart(xCons, yCons, n.reps = 1L, n.burn = c(2L, 1L, 1L)),
  pattern = "per-replication burn-in"
)

# and the consolidated names are gone from the signatures they left
for (name in c("resid.dist", "dispersion", "breaks", "max.rows")) {
  expect_false(name %in% names(formals(dbarts::bart)))
  expect_false(name %in% names(formals(dbarts::dbarts)))
}
expect_false("dart" %in% names(formals(dbarts::bart)))
expect_false("dart" %in% names(formals(dbarts::xbart)))
expect_false("levelGibbs" %in% names(formals(dbarts::bart)))

# the prior scalars dec-B116 moved onto the prior objects and the control
for (name in dbarts:::consolidatedPriorScalars) {
  expect_false(name %in% names(formals(dbarts::bart)))
}
expect_false("proposal.probs" %in% names(formals(dbarts::dbarts)))
# 'control' is a formal of both front doors again: xbart's refusal is reversed
expect_true("control" %in% names(formals(dbarts::bart)))
expect_true("control" %in% names(formals(dbarts::xbart)))
expect_true("proposal.probs" %in% names(formals(dbarts::dbartsControl)))

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
