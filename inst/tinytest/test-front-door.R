# The two front doors for the transition release: bart carries the modern
# interface and defaults, bartBT carries 0.9-34's argument list exactly, and
# a BayesTree-spelled bart call is forwarded to bartBT once per session.
# Also covers the renamed creation-time sigma estimate, the retired hurdle
# alias, and the refusal of a fit saved by 0.9-x.

resetWarnKey <- function(key) {
  env <- dbarts:::onceWarnState
  env[[key]] <- NULL
  invisible(NULL)
}

# every warning raised by an expression, by message, so an extra one cannot
# hide behind a pattern-only expectation
warningsOf <- function(expr) {
  seen <- character(0L)
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}

set.seed(2718)
nFD <- 90L
xFD <- matrix(runif(nFD * 3L), nFD, 3L)
colnames(xFD) <- c("x1", "x2", "x3")
linFD <- 2 * (xFD[, 1L] - 0.5) + xFD[, 2L] - xFD[, 3L]
yFD <- linFD + rnorm(nFD, 0, 0.5)
dfFD <- data.frame(xFD, y = yFD)

modernFD <- list(
  n.trees = 12L,
  n.samples = 40L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
legacyFD <- list(
  ntree = 12L,
  ndpost = 40L,
  nskip = 20L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)

# --- bartBT's signature is 0.9-34's, name for name and in order ---
expect_identical(
  names(formals(dbarts::bartBT)),
  c(
    "x.train",
    "y.train",
    "x.test",
    "sigest",
    "sigdf",
    "sigquant",
    "k",
    "power",
    "base",
    "splitprobs",
    "binaryOffset",
    "weights",
    "ntree",
    "ndpost",
    "nskip",
    "printevery",
    "keepevery",
    "keeptrainfits",
    "usequants",
    "numcut",
    "printcutoffs",
    "verbose",
    "nchain",
    "nthread",
    "combinechains",
    "keeptrees",
    "keepcall",
    "sampleronly",
    "seed",
    "proposalprobs",
    "keepsampler"
  )
)

# the names that select the legacy door: a derived set, pinned here so a
# formal added to either signature has to be looked at
expect_identical(
  sort(dbarts:::bartBTOnlyFormals),
  sort(c(
    "x.train",
    "y.train",
    "x.test",
    "splitprobs",
    "binaryOffset",
    "ntree",
    "ndpost",
    "nskip",
    "printevery",
    "keepevery",
    "keeptrainfits",
    "usequants",
    "numcut",
    "printcutoffs",
    "nchain",
    "nthread",
    "combinechains",
    "keeptrees",
    "keepcall",
    "sampleronly",
    "proposalprobs",
    "keepsampler"
  ))
)

# --- bart2 is a real formal list, not function(...) - a consumer reads its
# defaults - and forwards to bart bit for bit ---
expect_identical(formals(dbarts::bart2), formals(dbarts::bart))
expect_false(identical(names(formals(dbarts::bart2)), "..."))
expect_identical(
  eval(formals(dbarts::bart2)[["n.trees"]]),
  eval(formals(dbarts::bart)[["n.trees"]])
)

set.seed(41L)
fitModern <- do.call(dbarts::bart, c(list(y ~ x1 + x2 + x3, dfFD), modernFD))
resetWarnKey("tombstone.bart2")
set.seed(41L)
aliasWarnings <- character(0L)
fitAlias <- withCallingHandlers(
  do.call(dbarts::bart2, c(list(y ~ x1 + x2 + x3, dfFD), modernFD)),
  warning = function(w) {
    aliasWarnings <<- c(aliasWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_identical(fitAlias$yhat.train, fitModern$yhat.train)
expect_identical(fitAlias$sigma, fitModern$sigma)
expect_equal(length(aliasWarnings), 1L)
expect_true(grepl("'bart2' is now 'bart'", aliasWarnings[1L], fixed = TRUE))
# once per session: the second call is silent
expect_identical(
  length(warningsOf(
    do.call(dbarts::bart2, c(list(y ~ x1 + x2 + x3, dfFD), modernFD))
  )),
  0L
)

# --- a BayesTree-spelled bart call is the legacy door's fit, and warns
# exactly once per session ---
set.seed(31L)
fitLegacy <- do.call(dbarts::bartBT, c(list(xFD, yFD), legacyFD))
resetWarnKey("tombstone.bartShim")
set.seed(31L)
shimWarnings <- character(0L)
fitShim <- withCallingHandlers(
  do.call(dbarts::bart, c(list(x.train = xFD, y.train = yFD), legacyFD)),
  warning = function(w) {
    shimWarnings <<- c(shimWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_identical(fitShim$yhat.train, fitLegacy$yhat.train)
expect_identical(fitShim$sigma, fitLegacy$sigma)
expect_equal(length(shimWarnings), 1L)
expect_true(grepl("'bartBT'", shimWarnings[1L], fixed = TRUE))
expect_true(grepl("x.train", shimWarnings[1L], fixed = TRUE))
# the second BayesTree-spelled call in the session is silent
set.seed(31L)
expect_identical(
  length(warningsOf(
    fitShim2 <- do.call(dbarts::bart, c(list(xFD, yFD), legacyFD))
  )),
  0L
)
expect_identical(fitShim2$yhat.train, fitLegacy$yhat.train)

# the forwarded call is the one the caller wrote, so positional arguments
# stay positional and reach the legacy door's own names
expect_identical(dbarts:::callName(fitShim2$call), "bartBT")
expect_identical(
  dbarts:::callName(
    dbarts::bart(
      xFD,
      yFD,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 1L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE
    )$call
  ),
  "bart"
)

# --- a plain bart(x, y) with no BayesTree spelling is fit HERE, silently,
# at the MODERN defaults (75 trees, not 0.9-34's 200) ---
expect_identical(
  length(warningsOf(
    fitPlain <- dbarts::bart(xFD, yFD, samplerOnly = TRUE, verbose = FALSE)
  )),
  0L
)
expect_equal(fitPlain$control@n.trees, 75L)
expect_equal(
  dbarts::bartBT(xFD, yFD, sampleronly = TRUE, verbose = FALSE)$control@n.trees,
  200L
)

# --- the positional guard: 0.9-x's fourth positional argument was sigest,
# this one's is subset ---
expect_error(
  dbarts::bart(xFD, yFD, xFD[1:5, ], 1.0),
  pattern = "at most three positional arguments"
)
expect_error(
  dbarts::bart(xFD, yFD, xFD[1:5, ], 1.0),
  pattern = "bartBT"
)
# three positional arguments are the modern signature's own and are fine
expect_inherits(
  do.call(dbarts::bart, c(list(xFD, yFD, xFD[1:5, ]), modernFD)),
  "bart"
)

# --- '...' carries only registry names; anything else is a caller mistake ---
expect_error(
  do.call(dbarts::bart, c(list(xFD, yFD, zzzznotarg = 1), modernFD)),
  pattern = "unused argument 'zzzznotarg'"
)
expect_error(
  dbarts::dbartsControl(zzzznotarg = 1),
  pattern = "unused argument 'zzzznotarg'"
)
expect_error(
  dbarts::xbart(xFD, yFD, zzzznotarg = 1),
  pattern = "unused argument 'zzzznotarg'"
)

# rngSeed IS a registry name: accepted as seed, once-warned, and the value
# is used rather than dropped
resetWarnKey("tombstone.rngSeed.bart")
seedWarnings <- character(0L)
fitRngSeed <- withCallingHandlers(
  do.call(
    dbarts::bart,
    c(list(xFD, yFD, rngSeed = 77L, samplerOnly = TRUE), modernFD)
  ),
  warning = function(w) {
    seedWarnings <<- c(seedWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(seedWarnings), 1L)
expect_equal(fitRngSeed$control@seed, 77L)
resetWarnKey("tombstone.rngSeed.dbartsControl")
expect_equal(
  suppressWarnings(dbarts::dbartsControl(rngSeed = 88L))@seed,
  88L
)

# --- sigest is the creation-time estimate everywhere ---
expect_true("sigest" %in% names(formals(dbarts::dbarts)))
expect_true("sigest" %in% names(formals(dbarts::dbartsSpec)))
sigControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
expect_equal(
  dbarts::dbarts(xFD, yFD, control = sigControl, sigest = 1.25)$data@sigma,
  1.25
)
resetWarnKey("tombstone.sigma.dbarts")
sigmaWarnings <- character(0L)
sigmaSampler <- withCallingHandlers(
  dbarts::dbarts(xFD, yFD, control = sigControl, sigma = 1.25),
  warning = function(w) {
    sigmaWarnings <<- c(sigmaWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
expect_equal(length(sigmaWarnings), 1L)
expect_equal(sigmaSampler$data@sigma, 1.25)
expect_error(
  dbarts::dbarts(xFD, yFD, control = sigControl, sigma = 1.25, sigest = 2.0),
  pattern = "supply one"
)

# --- hurdle tokens: one spelling at the front door, none at dbarts() ---
expect_error(
  dbarts::dbarts(xFD, abs(yFD), family = "hurdle.lognormal"),
  pattern = "bart\\(x.train, y.train, family = \"hurdle.lognormal\"\\)"
)
twopartMsgFD <- tryCatch(
  dbarts::dbarts(xFD, abs(yFD), family = "twopart"),
  error = function(e) conditionMessage(e)
)
expect_true(grepl("hurdle.lognormal", twopartMsgFD, fixed = TRUE))
expect_error(
  do.call(
    dbarts::bart,
    c(list(xFD, abs(yFD), family = "twopart"), modernFD)
  ),
  pattern = "hurdle.lognormal"
)
expect_false("twopart" %in% eval(formals(dbarts::dbarts)[["family"]]))
expect_false(
  "hurdle.lognormal" %in% eval(formals(dbarts::dbarts)[["family"]])
)
expect_false("twopart" %in% eval(formals(dbarts::bart)[["family"]]))
expect_true(
  "hurdle.lognormal" %in% eval(formals(dbarts::bart)[["family"]])
)

# --- a fit saved by dbarts 0.9-x is refused by name ---
# built by hand: 0.9-x's state carried no format field at all, which is the
# whole version test. The fit object is otherwise a live one, so the refusal
# has to come from the state and not from a missing slot.
stateFit <- do.call(
  dbarts::bart,
  c(list(xFD, yFD, keepTrees = TRUE, keepSampler = TRUE), modernFD)
)
stateFit$fit$storeState()
legacyState <- stateFit$fit$state
attr(legacyState, "formatVersion") <- NULL
expect_error(
  stateFit$fit$setState(legacyState),
  pattern = "saved by dbarts 0.9-x"
)
# and through the re-creation branch getPointer takes after a load, which
# is what predict on a restored 0.9-x fit reaches: a round trip through
# serialization leaves the engine pointer dead, exactly as a save/load does
reloaded <- unserialize(serialize(stateFit, NULL))
reloaded$fit$state <- legacyState
expect_error(
  predict(reloaded, xFD[1:3, ]),
  pattern = "refit with this version"
)

rm(
  nFD,
  xFD,
  yFD,
  linFD,
  dfFD,
  modernFD,
  legacyFD,
  fitModern,
  fitAlias,
  fitLegacy,
  fitShim,
  fitShim2,
  fitPlain,
  fitRngSeed,
  sigControl,
  sigmaSampler,
  stateFit,
  reloaded,
  legacyState,
  aliasWarnings,
  shimWarnings,
  seedWarnings,
  sigmaWarnings,
  twopartMsgFD
)
