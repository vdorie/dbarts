# NULL is the spelling of "absent" on the public surface and NA is a missing
# value: an NA that 0.9-34 documented warns once per session until the
# tombstones expire, and one that never shipped is refused.

warnState <- dbarts:::onceWarnState
resetKeys <- function(...) {
  for (key in c(...)) {
    warnState[[key]] <- NULL
  }
}
naKey <- function(argument, caller) {
  paste0("tombstone.NA.", argument, ".", caller)
}
# the warnings an expression signals, muffled and returned
warningsOf <- function(expr) {
  seen <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}
expect_one_warning <- function(expr, pattern) {
  seen <- warningsOf(expr)
  expect_equal(length(seen), 1L)
  expect_true(grepl(pattern, seen[1L], fixed = TRUE))
  expect_true(grepl(dbarts:::tombstoneExpiry, seen[1L], fixed = TRUE))
}
expect_no_warning_of <- function(expr) {
  expect_equal(length(warningsOf(expr)), 0L)
}

set.seed(21L)
nObs <- 40L
x <- matrix(rnorm(nObs * 2L), nObs, 2L)
y <- x[, 1L] + rnorm(nObs)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 4L,
  n.burn = 2L,
  verbose = FALSE
)
bartArgs <- list(
  n.trees = 5L,
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
callBart <- function(...) {
  suppressMessages(do.call(dbarts::bart, c(list(x, y), bartArgs, list(...))))
}

# ---- site 1: sigest = NULL ----

for (formal in list(
  dbarts::bart,
  dbarts::dbarts,
  dbarts::dbartsSpec,
  dbarts::xbart
)) {
  expect_null(formals(formal)[["sigest"]])
}
expect_identical(formals(dbarts::bartBT)[["sigest"]], NA_real_)
expect_null(formals(dbarts::dbarts)[["sigma"]])
expect_null(formals(dbarts::dbartsSpec)[["sigma"]])

resetKeys(naKey("sigest", "bart"), naKey("sigest", "dbarts"))
expect_no_warning_of(callBart(sigest = NULL, seed = 3L))
defaultFit <- callBart(seed = 3L)
expect_identical(callBart(sigest = NULL, seed = 3L)$sigest, defaultFit$sigest)
expect_one_warning(
  naFit <- callBart(sigest = NA, seed = 3L),
  "'sigest = NA' is now 'sigest = NULL' on 'bart'"
)
# the NA is read as NULL, so the fit is the default one
expect_identical(naFit$sigest, defaultFit$sigest)
expect_identical(naFit$yhat.train, defaultFit$yhat.train)
expect_no_warning_of(callBart(sigest = NA, seed = 3L))
expect_error(callBart(sigest = -1), "'sigest'")
expect_error(callBart(sigest = "a"), "'sigest'")
expect_error(callBart(sigest = c(1, 2)), "'sigest'")
expect_identical(callBart(sigest = 2, seed = 3L)$sigest, 2)

sigmaSampler <- function(...) {
  dbarts::dbarts(x, y, control = control, ...)
}
defaultSigma <- sigmaSampler()$data@sigma
expect_no_warning_of(nullSampler <- sigmaSampler(sigest = NULL))
expect_identical(nullSampler$data@sigma, defaultSigma)
expect_one_warning(
  naSampler <- sigmaSampler(sigest = NA),
  "'sigest = NA' is now 'sigest = NULL' on 'dbarts'"
)
expect_identical(naSampler$data@sigma, defaultSigma)
expect_error(sigmaSampler(sigest = 0), "'sigest'")
expect_identical(sigmaSampler(sigest = 1.5)$data@sigma, 1.5)

# dbartsSpec: NULL leaves the data object's own estimate alone
specData <- dbarts::dbartsData(x, y)
specData@sigma <- 0.7
specOf <- function(...) {
  dbarts::dbartsSpec(specData, control = control, ...)$data@sigma
}
resetKeys(naKey("sigest", "dbartsSpec"))
expect_identical(specOf(), 0.7)
expect_identical(specOf(sigest = NULL), 0.7)
expect_one_warning(
  naSpec <- specOf(sigest = NA),
  "'sigest = NA' is now 'sigest = NULL' on 'dbartsSpec'"
)
expect_identical(naSpec, 0.7)
expect_identical(specOf(sigest = 1.5), 1.5)
expect_error(specOf(sigest = -1), "'sigest'")

xbartOf <- function(...) {
  dbarts::xbart(
    x,
    y,
    n.reps = 2L,
    n.trees = 5L,
    n.samples = 4L,
    n.burn = c(2L, 2L),
    n.test = 5,
    n.threads = 1L,
    seed = 5L,
    verbose = FALSE,
    ...
  )
}
resetKeys(naKey("sigest", "xbart"))
xbartDefault <- xbartOf()
expect_identical(xbartOf(sigest = NULL), xbartDefault)
expect_one_warning(
  naXbart <- xbartOf(sigest = NA),
  "'sigest = NA' is now 'sigest = NULL' on 'xbart'"
)
expect_identical(naXbart, xbartDefault)
expect_error(xbartOf(sigest = -1), "'sigest'")

# bartBT keeps BayesTree's NA silently
resetKeys(naKey("sigest", "bartBT"))
expect_no_warning_of(suppressMessages(dbarts::bartBT(
  x,
  y,
  sigest = NA,
  ndpost = 4L,
  nskip = 2L,
  ntree = 5L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)))

# the retired sigma spelling points at sigest and defaults to NULL
resetKeys("tombstone.sigma.dbarts", "tombstone.sigma.dbartsSpec")
resetKeys("tombstone.sigma.xbart", naKey("sigest", "dbarts"))
expect_no_warning_of(sigmaSampler(sigma = NULL))
expect_one_warning(
  viaSigma <- sigmaSampler(sigma = 1.5),
  "'sigma' is now 'sigest' on 'dbarts'"
)
expect_identical(viaSigma$data@sigma, 1.5)
# an NA under the old name is covered by the old name's own warning
resetKeys("tombstone.sigma.dbarts")
expect_one_warning(
  sigmaSampler(sigma = NA),
  "'sigma' is now 'sigest' on 'dbarts'"
)
expect_error(sigmaSampler(sigma = 1, sigest = 1), "supply one")
expect_one_warning(
  viaSigmaSpec <- specOf(sigma = 1.5),
  "'sigma' is now 'sigest' on 'dbartsSpec'"
)
expect_identical(viaSigmaSpec, 1.5)
expect_one_warning(
  viaSigmaXbart <- xbartOf(sigma = 1.5),
  "'sigma' is now 'sigest' on 'xbart'"
)
expect_identical(viaSigmaXbart, xbartOf(sigest = 1.5))
expect_error(xbartOf(sigma = 1, sigest = 1), "supply one")
sigmaEntries <- Filter(
  function(entry) identical(entry$name, "sigma"),
  dbarts:::dbartsTombstones
)
expect_true(all(vapply(
  sigmaEntries,
  function(entry) identical(entry$successor, "sigest"),
  NA
)))
expect_true(setequal(
  vapply(sigmaEntries, `[[`, "", "owner"),
  c("dbarts", "dbartsSpec", "xbart")
))
rm(defaultFit, naFit, nullSampler, naSampler, naSpec, xbartDefault, naXbart)
rm(viaSigma, viaSigmaSpec, viaSigmaXbart, sigmaEntries)

# ---- site 2: seed = NULL ----

seedCallers <- list(
  bart = function(...) callBart(...),
  bartBT = function(...) {
    suppressMessages(dbarts::bartBT(
      x,
      y,
      ndpost = 4L,
      nskip = 2L,
      ntree = 5L,
      nchain = 1L,
      nthread = 1L,
      verbose = FALSE,
      ...
    ))
  },
  dbarts = function(...) sigmaSampler(...),
  xbart = function(...) {
    dbarts::xbart(
      x,
      y,
      n.reps = 2L,
      n.trees = 5L,
      n.samples = 4L,
      n.burn = c(2L, 2L),
      n.test = 5,
      n.threads = 1L,
      verbose = FALSE,
      ...
    )
  },
  dbartsControl = function(...) dbarts::dbartsControl(...)
)
for (caller in names(seedCallers)) {
  expect_null(formals(get(caller, asNamespace("dbarts")))[["seed"]])
  resetKeys(naKey("seed", caller))
  expect_no_warning_of(seedCallers[[caller]](seed = NULL))
  expect_one_warning(
    seedCallers[[caller]](seed = NA),
    paste0("'seed = NA' is now 'seed = NULL' on '", caller, "'")
  )
  expect_no_warning_of(seedCallers[[caller]](seed = NA))
}
rm(caller)
# resolved to the value the default gives
expect_identical(
  dbarts::dbartsControl(seed = NULL)@seed,
  dbarts::dbartsControl()@seed
)
expect_identical(
  suppressWarnings(dbarts::dbartsControl(seed = NA)@seed),
  NA_integer_
)
seededControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 4L,
  seed = 9L
)
resetKeys(naKey("seed", "dbarts"))
expect_identical(
  suppressWarnings(
    dbarts::dbarts(x, y, control = seededControl, seed = NA)$control@seed
  ),
  9L
)
# the retired rngSeed = NA is covered by the retired name's own warning
resetKeys(naKey("seed", "dbartsControl"), "tombstone.rngSeed.dbartsControl")
expect_one_warning(
  dbarts::dbartsControl(rngSeed = NA),
  "'rngSeed' is now 'seed'"
)
# never shipped on these two: refused as a missing value, naming NULL
expect_error(
  dbarts::dbartsSpec(specData, control = control, seed = NA),
  "'seed' must not be NA on 'dbartsSpec'.*NULL"
)
expect_error(
  dbarts::dbartsValidateComposition(
    function() 0,
    function(theta) 0,
    function(theta, data) 0,
    function(state) state,
    function(state) 0,
    seed = NA
  ),
  "'seed' must not be NA on 'dbartsValidateComposition'.*NULL"
)
expect_error(
  dbarts::dbartsSpec(specData, control = control, seed = 1.5),
  "seed"
)

# ---- site 3: updateState = NULL ----

samplerMethods <- dbarts:::dbartsSampler$methods()
withUpdateState <- Filter(
  function(name) {
    method <- dbarts:::dbartsSampler$def@refMethods[[name]]
    is.function(method) && "updateState" %in% names(formals(method))
  },
  setdiff(samplerMethods, c("show", "callSuper", "copy", "initFields"))
)
expect_true(length(withUpdateState) >= 20L)
for (name in withUpdateState) {
  formal <- formals(dbarts:::dbartsSampler$def@refMethods[[name]])
  if (
    name %in%
      c("setTestPredictor", "setTestPredictorAndOffset", "setTestOffset")
  ) {
    next
  }
  expect_null(formal[["updateState"]], info = name)
}
expect_null(formals(dbarts::updatePredictorPerObservationJointly)$updateState)

resolve <- dbarts:::resolveUpdateState
onControl <- dbarts::dbartsControl(updateState = TRUE)
offControl <- dbarts::dbartsControl(updateState = FALSE)
expect_true(resolve(NULL, onControl))
expect_false(resolve(NULL, offControl))
expect_true(resolve(TRUE, offControl))
expect_false(resolve(FALSE, onControl))
resetKeys(naKey("updateState", "dbartsSampler"))
expect_one_warning(
  naResolved <- resolve(NA, offControl),
  "'updateState = NA' is now 'updateState = NULL' on 'dbartsSampler'"
)
# an NA reads as NULL, so it defers to the control
expect_false(naResolved)
expect_true(resolve(NA, onControl))
for (bad in list("yes", 1L, c(TRUE, FALSE), NA_character_[0L], list())) {
  expect_error(
    resolve(bad, onControl),
    "'updateState' must be TRUE, FALSE or NULL"
  )
}

mutable <- sigmaSampler()
sigmaBefore <- mutable$getSigmas()
expect_error(
  mutable$setSigma(3, updateState = "yes"),
  "'updateState' must be TRUE, FALSE or NULL"
)
# refused before the sampler changes
expect_identical(mutable$getSigmas(), sigmaBefore)
expect_no_warning_of(mutable$setSigma(3, updateState = NULL))
expect_equal(as.vector(mutable$getSigmas()), 3)
expect_error(
  mutable$run(0L, 1L, updateState = 1L),
  "'updateState' must be TRUE, FALSE or NULL"
)
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    mutable,
    x[, 1L],
    1L,
    updateState = "yes"
  ),
  "'updateState' must be TRUE, FALSE or NULL"
)
resetKeys(naKey("updateState", "dbartsSampler"))
expect_one_warning(
  mutable$setSigma(2, updateState = NA),
  "'updateState = NA' is now 'updateState = NULL' on 'dbartsSampler'"
)
rm(mutable, sigmaBefore, naResolved)

# ---- site 4: treeShift ----

expect_identical(
  eval(formals(dbarts::dbartsControl)[["treeShift"]]),
  c("auto", "always", "never")
)
expect_false("levelGibbs" %in% names(formals(dbarts::dbartsControl)))
expect_false("levelGibbs" %in% names(formals(dbarts:::cgm)))
expect_false("levelGibbs" %in% names(formals(dbarts:::dart)))
# the bridge keeps reading the tri-state logical slot, which round-trips
for (word in c("auto", "always", "never")) {
  expect_identical(
    dbarts:::controlArgumentFromSlot(
      "treeShift",
      dbarts::dbartsControl(treeShift = word)
    ),
    word
  )
}
expect_identical(
  dbarts::dbartsControl(treeShift = "always")@levelGibbs,
  TRUE
)
expect_identical(dbarts::dbartsControl()@levelGibbs, NA)
# a control's treeShift survives a front door that rebuilds the control
fitTreeShift <- callBart(
  control = dbarts::dbartsControl(treeShift = "never"),
  keepSampler = TRUE
)
expect_false(fitTreeShift$fit$control@levelGibbs)
rm(fitTreeShift)

# ---- site 5: nbinom(dispersion = NULL) ----

nbinom <- dbarts:::dbartsFamilies$nbinom
expect_null(formals(nbinom)[["dispersion"]])
expect_true(is.na(nbinom()@settings$dispersion))
expect_identical(nbinom(dispersion = NULL)@settings, nbinom()@settings)
expect_equal(nbinom(dispersion = 3)@settings$dispersion, 3)
for (bad in list(NA, NA_real_, 0, -1, Inf, "a", c(1, 2))) {
  expect_error(nbinom(dispersion = bad), "'dispersion' must be NULL")
}

# ---- site 6: run(numBurnIn = NULL, numSamples = NULL) ----

runner <- sigmaSampler()
expect_null(formals(dbarts:::dbartsSampler$def@refMethods[["run"]])$numBurnIn)
expect_null(formals(dbarts:::dbartsSampler$def@refMethods[["run"]])$numSamples)
byControl <- runner$run()
expect_equal(dim(byControl$train), c(nObs, 4L))
expect_identical(dim(runner$run(NULL, NULL)$train), dim(byControl$train))
expect_identical(dim(runner$run(0L, NULL)$train), dim(byControl$train))
resetKeys(
  naKey("numBurnIn", "dbartsSampler"),
  naKey("numSamples", "dbartsSampler")
)
seen <- warningsOf(naRun <- runner$run(NA, NA))
expect_equal(length(seen), 2L)
expect_true(any(grepl(
  "'numBurnIn = NA' is now 'numBurnIn = NULL'",
  seen,
  fixed = TRUE
)))
expect_true(any(grepl(
  "'numSamples = NA' is now 'numSamples = NULL'",
  seen,
  fixed = TRUE
)))
expect_identical(dim(naRun$train), dim(byControl$train))
expect_equal(dim(runner$run(0L, 2L)$train), c(nObs, 2L))
expect_error(runner$run("a", 2L), "numBurnIn")
expect_error(runner$run(0L, 1.5), "numSamples")
rm(runner, byControl, naRun, seen)

# ---- site 7: dbartsControl(n.samples = NULL) ----

expect_null(formals(dbarts::dbartsControl)[["n.samples"]])
expect_identical(dbarts::dbartsControl()@n.samples, NA_integer_)
expect_identical(dbarts::dbartsControl(n.samples = NULL)@n.samples, NA_integer_)
resetKeys(naKey("n.samples", "dbartsControl"))
expect_one_warning(
  naControl <- dbarts::dbartsControl(n.samples = NA),
  "'n.samples = NA' is now 'n.samples = NULL' on 'dbartsControl'"
)
expect_identical(naControl@n.samples, NA_integer_)
expect_identical(dbarts::dbartsControl(n.samples = 7L)@n.samples, 7L)
expect_error(dbarts::dbartsControl(n.samples = "a"), "n.samples")
expect_error(dbarts::dbartsControl(n.samples = -1L), "n.samples")
rm(naControl)

# ---- site 8: NA refused where NULL is the default ----

priors <- dbarts:::dbartsPriors
for (argument in c("sd", "scale")) {
  for (constructor in c("normal", "linear", "gp")) {
    args <- if (constructor == "normal") list() else list(columns = 1L)
    args[[argument]] <- NA
    expect_error(
      do.call(priors[[constructor]], args),
      paste0("'", argument, "' must not be NA"),
      info = paste(constructor, argument)
    )
    args[[argument]] <- NULL
    expect_silent(do.call(priors[[constructor]], args))
  }
}
expect_identical(priors$normal(sd = NULL)@prior.sd, NA_real_)
expect_error(priors$normal(sd = -1), "'sd' must be positive")
expect_error(priors$normal(sd = "a"), "'sd'")
expect_error(priors$dart(rho = NA), "'rho' must not be NA")
expect_error(priors$dart(update.delay = NA), "'update.delay' must not be NA")
expect_error(priors$dart(rho = -1), "rho")
expect_error(priors$dart(update.delay = -1), "update.delay")
expect_identical(priors$dart(rho = NULL)@rho, NA_real_)
expect_identical(priors$dart(update.delay = NULL)@update.delay, NA_real_)

# ---- site 9: the 0.9-34 test setters keep an updateState ----

setterSampler <- dbarts::dbarts(x, y, test = x[1:5, ], control = control)
resetKeys("tombstone.testSetters.updateState")
expect_no_warning_of(setterSampler$setTestPredictor(x[1:5, ]))
expect_no_warning_of(setterSampler$setTestOffset(NULL, NULL))
expect_one_warning(
  setterSampler$setTestPredictor(x[1:5, ], updateState = TRUE),
  "the test-data setters no longer take 'updateState'"
)
# once per session, whichever setter carries it
expect_no_warning_of(setterSampler$setTestOffset(NULL, updateState = FALSE))
expect_no_warning_of(
  setterSampler$setTestPredictorAndOffset(x[1:5, ], NULL, updateState = NA)
)
expect_no_warning_of(setterSampler$setTestPredictorAndOffset(x[1:5, ], NULL))
expect_equal(dim(setterSampler$data@x.test), c(5L, 2L))
rm(setterSampler)
registered <- vapply(dbarts:::dbartsTombstones, `[[`, "", "name")
expect_true("NA sigest" %in% registered)
expect_true("NA seed" %in% registered)
expect_true("updateState on the test setters" %in% registered)

# ---- site 10: gaussian(link = identity) unquoted ----

gaussianFamily <- dbarts:::dbartsFamilies$gaussian
expect_identical(gaussianFamily(link = identity), gaussianFamily())
expect_identical(gaussianFamily(link = "identity"), gaussianFamily())
identityName <- "identity"
expect_identical(gaussianFamily(identityName), gaussianFamily())
expect_error(gaussianFamily(link = log), "not supported")
expect_error(gaussianFamily(link = "log"), "not supported")
expect_error(gaussianFamily(link = 3), "takes 'link' first")
gaussianFit <- callBart(family = gaussian(link = identity), seed = 3L)
expect_identical(
  gaussianFit$yhat.train,
  callBart(family = "gaussian", seed = 3L)$yhat.train
)
