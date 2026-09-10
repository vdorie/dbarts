# dec-B116: bart() sheds the seven scalars that mirror slots on the prior
# objects and the tree-move mixture, each riding '...' for one release, and
# both front doors gain 'control ='. Covers the retirements (warned once,
# draws identical to the object spelling), the collision refusals, the four
# engine settings reached through a control from either door, the
# flat-wins-over-the-slot rule, and the fit-state refusal.

resetWarn <- function(key) {
  env <- dbarts:::onceWarnState
  env[[key]] <- NULL
  invisible(NULL)
}

# every warning an expression raises, by message: an extra one cannot hide
# behind a pattern-only expectation
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

set.seed(3141L)
nFC <- 60L
xFC <- matrix(runif(nFC * 3L), nFC, 3L, dimnames = list(NULL, c("a", "b", "c")))
yFC <- xFC[, 1L] - xFC[, 3L] + rnorm(nFC, 0, 0.3)

quickFC <- list(
  n.trees = 5L,
  n.samples = 10L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 17L
)
fitFC <- function(...) {
  do.call(dbarts::bart, c(list(xFC, yFC), quickFC, list(...)))
}
drawsOf <- function(fit) list(fit$yhat.train, fit$sigma)

# ---- each retired name warns once and fits what the object spelling fits ---

# name, the retired call's arguments, and the object spelling that replaces it
retirements <- list(
  power = list(
    old = list(power = 3.0),
    new = list(tree.prior = quote(dbarts::dbartsPriors$cgm(3.0, 0.95)))
  ),
  base = list(
    old = list(base = 0.8),
    new = list(tree.prior = quote(dbarts::dbartsPriors$cgm(2.0, 0.8)))
  ),
  split.probs = list(
    old = list(split.probs = c(0.5, 0.25, 0.25)),
    new = list(
      tree.prior = quote(dbarts::dbartsPriors$cgm(
        2.0,
        0.95,
        c(0.5, 0.25, 0.25)
      ))
    )
  ),
  prior.scale = list(
    old = list(prior.scale = 1.5),
    new = list(node.prior = quote(dbarts::dbartsPriors$normal(scale = 1.5)))
  ),
  sigdf = list(
    old = list(sigdf = 5.0),
    new = list(resid.prior = quote(dbarts::dbartsPriors$chisq(5.0, 0.90)))
  ),
  sigquant = list(
    old = list(sigquant = 0.75),
    new = list(resid.prior = quote(dbarts::dbartsPriors$chisq(3.0, 0.75)))
  ),
  proposal.probs = list(
    old = list(proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4)),
    new = list(
      control = quote(dbarts::dbartsControl(
        proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4)
      ))
    )
  )
)

for (name in names(retirements)) {
  resetWarn(paste0("tombstone.consolidated.", name, ".bart"))
  spelling <- retirements[[name]]
  retiredDraws <- NULL
  warned <- warningsOf(
    retiredDraws <- drawsOf(do.call(fitFC, spelling$old))
  )
  expect_equal(length(warned), 1L, info = name)
  expect_true(
    grepl(paste0("'", name, "' has left 'bart'"), warned[1L]),
    info = name
  )
  # once per session: the second call is silent and still applies the value
  expect_equal(
    length(warningsOf(secondDraws <- drawsOf(do.call(fitFC, spelling$old)))),
    0L,
    info = name
  )
  expect_identical(retiredDraws, secondDraws, info = name)
  objectDraws <- drawsOf(do.call(
    fitFC,
    lapply(spelling$new, function(e) if (is.call(e)) eval(e) else e)
  ))
  expect_identical(retiredDraws, objectDraws, info = name)
}
rm(retirements, retiredDraws, secondDraws, objectDraws, warned, name, spelling)

# ---- the collision refusals ------------------------------------------------

expect_error(
  suppressWarnings(fitFC(tree.prior = dbarts::dbartsPriors$cgm(), power = 3.0)),
  "cannot be combined with 'power'"
)
expect_error(
  suppressWarnings(fitFC(tree.prior = dbarts::dbartsPriors$cgm(), base = 0.8)),
  "cannot be combined with 'base'"
)
expect_error(
  suppressWarnings(fitFC(
    tree.prior = dbarts::dbartsPriors$cgm(),
    split.probs = c(0.5, 0.25, 0.25)
  )),
  "cannot be combined with 'split.probs'"
)
expect_error(
  suppressWarnings(fitFC(
    node.prior = dbarts::dbartsPriors$normal(2),
    prior.scale = 1.5
  )),
  "cannot be combined with 'prior.scale'"
)
expect_error(
  suppressWarnings(fitFC(
    resid.prior = dbarts::dbartsPriors$chisq(),
    sigdf = 5.0
  )),
  "cannot be combined with 'sigdf'"
)
expect_error(
  suppressWarnings(fitFC(
    resid.prior = dbarts::dbartsPriors$chisq(),
    sigquant = 0.75
  )),
  "cannot be combined with 'sigquant'"
)

# ---- control reaches the four engine settings from both doors --------------

engineControl <- dbarts::dbartsControl(
  categoricalExhaustiveCap = 7L,
  testFitParallelCutoff = 111L,
  predictParallelCutoff = 222L,
  sparseDensityThreshold = 0.35
)
engineFit <- fitFC(control = engineControl, keepSampler = TRUE)
expect_equal(engineFit$fit$control@categoricalExhaustiveCap, 7L)
expect_equal(engineFit$fit$control@testFitParallelCutoff, 111L)
expect_equal(engineFit$fit$control@predictParallelCutoff, 222L)
expect_equal(engineFit$fit$control@sparseDensityThreshold, 0.35)
# and every setting bart carries a flat name for is still bart's own default,
# not the control's
expect_equal(engineFit$fit$control@n.trees, 5L)
expect_equal(engineFit$fit$control@n.chains, 1L)

# the same control on xbart's own sweep, read off the control its cells are
# built from
# the tracer runs in the traced function's own frame, so the capture is an
# environment handed to it rather than a name reached by '<<-'
captured <- new.env(parent = emptyenv())
trace(
  dbarts:::bartcoreDataHandle,
  tracer = bquote(assign("control", control, envir = .(captured))),
  print = FALSE
)
xval <- dbarts::xbart(
  xFC,
  yFC,
  n.samples = 5L,
  n.burn = c(3L, 2L),
  n.reps = 1L,
  n.trees = 5L,
  k = c(1, 2),
  n.threads = 1L,
  seed = 19L,
  control = engineControl
)
untrace(dbarts:::bartcoreDataHandle)
xbartControl <- captured$control
# xbart's grid still runs
expect_equal(dim(xval), c(1L, 2L))
expect_true(all(is.finite(xval)))
expect_equal(xbartControl@categoricalExhaustiveCap, 7L)
expect_equal(xbartControl@testFitParallelCutoff, 111L)
expect_equal(xbartControl@predictParallelCutoff, 222L)
expect_equal(xbartControl@sparseDensityThreshold, 0.35)
# the fields the sweep forces stand whatever the control carried
expect_equal(xbartControl@n.chains, 1L)
expect_equal(xbartControl@n.threads, 1L)
expect_false(xbartControl@keepTrees)
expect_false(xbartControl@keepTrainingFits)
expect_false(xbartControl@updateState)
expect_false(xbartControl@verbose)
rm(xval, xbartControl, captured)

# ---- flat wins over the slot; an unnamed one takes the slot ----------------

flatControl <- dbarts::dbartsControl(n.cuts = 40L, n.thin = 2L)
flatFit <- fitFC(control = flatControl, n.cuts = 25L, keepSampler = TRUE)
expect_equal(flatFit$fit$control@n.cuts, 25L)
expect_equal(flatFit$fit$control@n.thin, 2L)

# ---- a control carrying fit state is refused by name -----------------------

varianceSampler <- dbarts::dbarts(
  xFC,
  yFC,
  variance = ~a,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 5L,
    updateState = FALSE
  )
)
expect_true(any(startsWith(
  names(attributes(varianceSampler$control)),
  "bartcore."
)))
expect_error(
  fitFC(control = varianceSampler$control),
  "pass a fresh dbartsControl"
)
expect_error(
  dbarts::xbart(
    xFC,
    yFC,
    n.reps = 1L,
    n.threads = 1L,
    control = varianceSampler$control
  ),
  "pass a fresh dbartsControl"
)
rm(varianceSampler)

# ---- the mixture is still settable between runs, now through setControl ----

midRun <- dbarts::dbarts(
  xFC,
  yFC,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 3L,
    n.burn = 0L,
    n.samples = 5L,
    updateState = FALSE
  )
)
invisible(midRun$run())
frozenNow <- midRun$control
frozenNow@proposal.probs[
  c("birth_death", "swap", "change", "perturb", "rule_gibbs")
] <- 0
midRun$setControl(frozenNow)
expect_equal(midRun$control@proposal.probs[["birth_death"]], 0)
expect_true(all(is.finite(midRun$run()$train)))
rm(midRun, frozenNow)

# ---- the mixture's validity travels with the slot ---------------------------

expect_error(
  dbarts::dbartsControl(
    proposal.probs = c(birth_death = 0.7, swap = 0, change = 0.4)
  ),
  "sum to 1"
)
expect_error(
  dbarts::dbartsControl(proposal.probs = c(birth_death = 0.7, change = 0.4)),
  "must be in \\[0, 1\\]"
)
expect_error(
  dbarts::dbartsControl(proposal.probs = c(swap = 0.1)),
  "name at least one of"
)
expect_error(
  dbarts::dbartsControl(proposal.probs = c(birth_death = 1, birth = 0)),
  "birth probability"
)
# a partial spelling is filled to the canonical six, in order
expect_identical(
  names(
    dbarts::dbartsControl(proposal.probs = c(birth_death = 0.7))@proposal.probs
  ),
  c("birth_death", "swap", "change", "perturb", "rule_gibbs", "birth")
)
# and a slot rewritten past validity is refused when the sampler is built
badControl <- dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
badControl@proposal.probs <- badControl@proposal.probs[1:3]
expect_error(dbarts::dbarts(xFC, yFC, control = badControl), "proposal.probs")
rm(badControl)

# ---- a fit saved with the old model slots still loads -----------------------

# the mixture used to live on dbartsModel; a saved object still carries those
# attributes, which the class no longer declares. They must be inert, not a
# validity failure or a refused install.
legacySampler <- dbarts::dbarts(
  xFC,
  yFC,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 3L,
    n.samples = 3L,
    updateState = FALSE
  )
)
legacyModel <- legacySampler$model
attr(legacyModel, "p.birth_death") <- 1.0
attr(legacyModel, "p.swap") <- 0.0
attr(legacyModel, "p.change") <- 0.0
attr(legacyModel, "p.perturb") <- 0.0
attr(legacyModel, "p.rule_gibbs") <- 0.0
attr(legacyModel, "p.birth") <- 0.5
expect_true(isTRUE(validObject(legacyModel, test = TRUE)))
expect_silent(legacySampler$setModel(legacyModel))
# and the mixture in force is still the control's, not the stale attributes'
expect_equal(legacySampler$control@proposal.probs[["birth_death"]], 0.6)
expect_true(all(is.finite(legacySampler$run()$train)))
rm(legacySampler, legacyModel)
