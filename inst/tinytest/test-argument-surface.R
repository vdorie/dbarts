# The bart2 argument-consolidation arc's maintained contracts
# (docs/plans/bart2-argument-consolidation.md 4.6): T-D (S2, family gating),
# T-B (S3, shared-default text), and T-A (S5, dbartsControl formal parity)
# already live here; T-C (bart -> bart2 concept map) and T-E (the per-forest
# reconstruction identity) arrive with S10, S11 respectively.

# T-D. One diagnosis test per row of 3.c.2: the classed warning fires under
# a family that ignores the argument and not under one that uses it; plus
# monotonicity assertions (N7) that already-loud refusals keep their
# severity and message.

set.seed(101)
n <- 40L
x <- matrix(rnorm(n * 2L), n, dimnames = list(NULL, c("a", "b")))
quick <- list(
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

y.gaussian <- rnorm(n)
y.binary <- rbinom(n, 1L, 0.5)
y.multi <- factor(sample(letters[1:3], n, replace = TRUE))
y.ordinal <- factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
y.count <- rpois(n, 4)
y.aft <- cbind(abs(rnorm(n)) + 0.1, rep(1L, n))

fit2 <- function(y, ...) {
  do.call(dbarts::bart2, c(list(x, y), quick, list(...)))
}

# one inventory row x one representative family/site each, class asserted
expect_warning(
  fit2(y.multi, family = "multinomial", sigest = 5),
  pattern = "sigest",
  class = "dbartsFamilyGatedWarning"
)
expect_warning(
  fit2(y.ordinal, family = "ordinal", dispersion = 2),
  pattern = "dispersion",
  class = "dbartsFamilyGatedWarning"
)
expect_warning(
  fit2(y.count, family = "nbinom", breaks = 5),
  pattern = "breaks",
  class = "dbartsFamilyGatedWarning"
)

# counter-tests: live families stay silent
expect_silent(fit2(y.gaussian, family = "gaussian", sigest = 5))
expect_silent(fit2(y.aft, family = "aft", sigest = 5))

# the suppliedness bit: a defaulted (unsupplied) argument never gates
expect_silent(fit2(y.binary, family = "probit"))

# two gated names supplied together still emit exactly one warning
countWarnings <- function(expr, class) {
  count <- 0L
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, class)) {
        count <<- count + 1L
      }
      invokeRestart("muffleWarning")
    }
  )
  count
}
expect_equal(
  countWarnings(
    fit2(y.binary, family = "probit", sigest = 5, dispersion = 2),
    "dbartsFamilyGatedWarning"
  ),
  1L
)

# samplerOnly = TRUE still warns (the standard-site reordering)
expect_warning(
  fit2(y.binary, family = "probit", sigest = 5, samplerOnly = TRUE),
  class = "dbartsFamilyGatedWarning"
)

# hurdle.lognormal: sigest is live on the positive half, so the rule's own
# scoping ("an argument whose only effect is on a family this fit is not")
# makes a diagnosis against the occupancy half's forced "probit" a FALSE
# one; bart2Hurdle strips it from that component call rather than let it
# leak. No warning at all - not even the occupancy component's.
y.hurdle <- c(rep(0, n / 2L), abs(rnorm(n / 2L)) + 0.1)
expect_equal(
  countWarnings(
    fit2(y.hurdle, family = "hurdle.lognormal", sigest = 5),
    "dbartsFamilyGatedWarning"
  ),
  0L
)

# dispersion is inert on both components; the outer site diagnoses it once
# and bart2Hurdle strips it from both component calls, so exactly one
# warning reaches the caller, naming the outer family
expect_equal(
  countWarnings(
    fit2(y.hurdle, family = "hurdle.lognormal", dispersion = 2),
    "dbartsFamilyGatedWarning"
  ),
  1L
)
expect_warning(
  fit2(y.hurdle, family = "hurdle.lognormal", dispersion = 2),
  pattern = "hurdle.lognormal",
  class = "dbartsFamilyGatedWarning"
)

# bart() never calls the helper - silence is preserved by construction
expect_silent(dbarts::bart(
  x,
  y.binary,
  sigest = 5,
  ntree = 3L,
  ndpost = 5L,
  nskip = 2L,
  verbose = FALSE
))

# rbart_vi warns on its resolved probit path
group <- factor(rep(1:4, length.out = n))
expect_warning(
  dbarts::rbart_vi(
    x,
    factor(y.binary),
    group.by = group,
    sigest = 5,
    n.thin = 1L,
    n.trees = 3L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  class = "dbartsFamilyGatedWarning"
)

# monotonicity (N7): already-loud refusals keep their severity and message
expect_error(
  fit2(y.multi, family = "multinomial", samplerOnly = TRUE),
  pattern = "does not support 'samplerOnly'"
)
expect_error(
  fit2(y.multi, family = "multinomial", weights = rep(1, n)),
  pattern = "does not support 'weights'"
)
expect_error(
  fit2(y.multi, family = "multinomial", prior.scale = 2),
  pattern = "prior.scale"
)

# T-B (S3). For every name shared by bart2 and dbarts, the deparsed default
# expressions agree, except the table below - copied from the plan (4.6).
# tree.prior/node.prior/resid.prior became bart2 formals at S7, so the loop
# now walks all three; they stay excepted because bart2's NULL means "build
# from the shorthands" while dbarts's own defaults are the bare constructors.
tbExceptions <- data.frame(
  name = c(
    "verbose",
    "n.samples",
    "family",
    "tree.prior",
    "node.prior",
    "resid.prior"
  ),
  reason = c(
    "fitters announce, constructors do not",
    "different roles; semantics differ (d2)",
    "multinomial is a bart2-only composition",
    "bart2's NULL means \"build from the shorthands\"",
    "same",
    "same"
  ),
  stringsAsFactors = FALSE
)
bart2Formals <- formals(dbarts::bart2)
dbartsFormals <- formals(dbarts::dbarts)
sharedFormalNames <- setdiff(
  intersect(names(bart2Formals), names(dbartsFormals)),
  tbExceptions$name
)
diverged <- Filter(
  function(nm) {
    !identical(deparse(bart2Formals[[nm]]), deparse(dbartsFormals[[nm]]))
  },
  sharedFormalNames
)
expect_equal(diverged, character(0))

# match.arg error messages for bad tokens: bart2 resolves factors/missing in
# its own frame before forwarding (S3), so a bad token errors here, naming
# the choices, same as R's own match.arg does for any other formal
expect_error(fit2(y.gaussian, factors = "bogus"), pattern = "should be one of")
expect_error(fit2(y.gaussian, missing = "bogus"), pattern = "should be one of")

# S2 interaction: the suppliedNames snapshot (argNames, R/bart.R) precedes
# the S3 resolve-and-forward step, so making factors/missing/proposal.probs
# unconditionally forwarded does not make them look "supplied" to family
# gating - a defaulted factors/missing call under a gating family stays
# silent
expect_silent(fit2(y.multi, family = "multinomial"))
expect_silent(fit2(y.count, family = "nbinom"))

# defaulted-vs-explicit equivalence (S3, D6): the resolved default is
# exactly what dbarts() already applied, so forwarding it explicitly changes
# no draw
sameDraws <- function(a, b) {
  identical(a$yhat.train, b$yhat.train) && identical(a$sigma, b$sigma)
}
defaultedProbs <- fit2(y.gaussian, seed = 77L)
explicitProbs <- fit2(
  y.gaussian,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  seed = 77L
)
expect_true(sameDraws(defaultedProbs, explicitProbs))

dfFactor <- data.frame(
  x1 = rnorm(n),
  x2 = factor(sample(letters[1:3], n, replace = TRUE))
)
defaultedFactors <- do.call(
  dbarts::bart2,
  c(list(y.gaussian ~ ., dfFactor), quick, list(seed = 77L))
)
explicitFactors <- do.call(
  dbarts::bart2,
  c(
    list(y.gaussian ~ ., dfFactor),
    quick,
    list(factors = "categorical", seed = 77L)
  )
)
expect_true(sameDraws(defaultedFactors, explicitFactors))

defaultedMissing <- fit2(y.gaussian, seed = 77L)
explicitMissing <- fit2(y.gaussian, missing = "incorporate", seed = 77L)
expect_true(sameDraws(defaultedMissing, explicitMissing))

# D4's second half: proposal.probs' default is now the named vector, always
# forwarded - a defaulted bart2 call still composes with monotone
expect_silent(fit2(y.gaussian, monotone = c(a = "+")))

# rbart_vi shares the pattern for factors/missing only (no proposal.probs
# formal)
defaultedRbart <- dbarts::rbart_vi(
  x,
  y.gaussian,
  group.by = group,
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  n.thin = 1L,
  verbose = FALSE,
  seed = 77L
)
explicitRbart <- dbarts::rbart_vi(
  x,
  y.gaussian,
  group.by = group,
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  n.thin = 1L,
  verbose = FALSE,
  seed = 77L,
  factors = "categorical",
  missing = "incorporate"
)
expect_true(sameDraws(defaultedRbart, explicitRbart))

# S4 (docs/plans/bart2-argument-consolidation.md 3.e, 4.1, 4.4): storage/
# updateState are now real formals, appended at the end of both entry
# points; '...' is rejection-only, diagnosed by a shared helper.

lastTwoNamed <- function(fn) {
  fnFormals <- names(formals(fn))
  fnFormals[length(fnFormals) - c(2L, 1L)]
}
expect_equal(lastTwoNamed(dbarts::bart2), c("storage", "updateState"))
expect_equal(lastTwoNamed(dbarts::rbart_vi), c("storage", "updateState"))

# storage/updateState reach the control with explicit, non-default values
explicitControl <- fit2(
  y.gaussian,
  storage = "single",
  updateState = FALSE,
  samplerOnly = TRUE
)
expect_equal(explicitControl$control@storage, "single")
expect_equal(explicitControl$control@updateState, FALSE)
defaultControl <- fit2(y.gaussian, samplerOnly = TRUE)
expect_equal(defaultControl$control@storage, "double")
expect_equal(defaultControl$control@updateState, TRUE)

rbartControlFit <- dbarts::rbart_vi(
  x,
  y.gaussian,
  group.by = group,
  storage = "single",
  updateState = FALSE,
  seed = 314L,
  keepSampler = TRUE,
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  n.thin = 1L,
  verbose = FALSE
)
rbartControl <- rbartControlFit$fit[[1L]]$control
expect_equal(rbartControl@storage, "single")
expect_equal(rbartControl@updateState, FALSE)
expect_equal(rbartControl@seed, 314L)

# dots rejection: an unknown name errors naming itself and suggesting the
# nearest formal; a name with no close match still errors, unsuggested
expect_error(fit2(y.gaussian, n.tres = 5), pattern = "did you mean 'n.trees'")
expect_error(
  dbarts::rbart_vi(x, y.gaussian, group.by = group, n.tres = 5),
  pattern = "did you mean 'n.trees'"
)
expect_error(fit2(y.gaussian, zzzznotarg = 5), pattern = "zzzznotarg")

# rngSeed (S5, b1/fork 5): dbartsControl's own spelling through S4; the
# passthrough that let it keep flowing through '...' is gone now that the
# rename lands, and it is the first retiredDotsNames row instead - a named
# refusal pointing at the rename, not a nearest-formal guess
expect_error(
  fit2(y.gaussian, rngSeed = 99L, samplerOnly = TRUE),
  pattern = "unknown argument 'rngSeed': 'rngSeed' was renamed to 'seed'"
)
expect_error(
  dbarts::rbart_vi(x, y.gaussian, group.by = group, rngSeed = 99L),
  pattern = "unknown argument 'rngSeed': 'rngSeed' was renamed to 'seed'"
)

# partial matching still works for a real formal
partialFit <- do.call(
  dbarts::bart2,
  c(
    list(x, y.gaussian, n.sampl = 5L, samplerOnly = TRUE),
    quick[setdiff(names(quick), "n.samples")]
  )
)
expect_equal(partialFit$control@n.samples, 5L)

# T-A (S5, 4.6): every dbartsControl formal is a bart2 formal, spelled
# identically - an explicit list of the 16 slots as of 1.0-0, a freeze-time
# snapshot rather than a ratchet (10.7 r2 M8): a dbartsControl formal added
# after 1.0-0 does not silently enlarge what this test checks.
controlFormals1_0_0 <- c(
  "verbose",
  "keepTrainingFits",
  "useQuantiles",
  "keepTrees",
  "storage",
  "n.samples",
  "n.cuts",
  "n.burn",
  "n.trees",
  "n.chains",
  "n.threads",
  "n.thin",
  "printEvery",
  "printCutoffs",
  "seed",
  "updateState"
)
expect_equal(
  sort(names(formals(dbarts::dbartsControl))),
  sort(controlFormals1_0_0)
)
expect_true(setequal(
  intersect(controlFormals1_0_0, names(formals(dbarts::bart2))),
  controlFormals1_0_0
))
expect_true(setequal(
  intersect(controlFormals1_0_0, names(formals(dbarts::rbart_vi))),
  controlFormals1_0_0
))

# S6 (docs/plans/bart2-argument-consolidation.md 3.b b3, 4.2, 6.3): the
# variance quartet collapses to a dedicated varianceForest() constructor;
# variance = keeps its shorthand (NULL/FALSE/TRUE/formula/character/index)
# and additionally accepts a varianceForest object (fork 2b). The flat
# n.trees.variance/power.variance/base.variance formals are gone from
# bart2/dbarts/dbartsSpec; the base-vs-HEAD byte-identity gate for the
# removed flat spelling rides an out-of-tree A/B probe (a second installed
# build is needed to exercise a formal this build no longer has), recorded
# in the landing notes. What lives here runs entirely against HEAD.

varN <- 40L
varX <- matrix(
  rnorm(varN * 2L),
  varN,
  dimnames = list(NULL, c("x1", "x2"))
)
varY <- rnorm(varN)
varControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 2L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE
)
varianceAttr <- function(variance, factors = "categorical") {
  attr(
    dbarts::dbarts(
      varX,
      varY,
      control = varControl,
      variance = variance,
      factors = factors
    )$control,
    "bartcore.variance"
  )
}

# round-trip identity (b3, one argument two accepted types): a varianceForest
# with no knobs declared resolves BYTE-IDENTICALLY to the shorthand it wraps,
# both when it reads every column and when it is restricted
expect_identical(varianceAttr(TRUE), varianceAttr(dbarts::varianceForest()))
expect_identical(
  varianceAttr(~x1),
  varianceAttr(dbarts::varianceForest(vars = ~x1))
)

# a declared n.trees/base/power on the object is honored and distinguishes it
# from the plain shorthand's own defaults
objectAttr <- varianceAttr(
  dbarts::varianceForest(vars = ~x1, n.trees = 20L, base = 0.9, power = 1.5)
)
expect_equal(objectAttr$n.trees, 20L)
expect_equal(objectAttr$base, 0.9)
expect_equal(objectAttr$power, 1.5)
expect_false(identical(objectAttr, varianceAttr(~x1)))

# factor-term expansion through vars = ~z under factors = "indicators" (the
# b3 selector requirement - vars resolves through the SAME resolveVarianceColumns
# the shorthand uses, so a factor term expands identically through either route)
dfFactorVar <- data.frame(
  x1 = rnorm(varN),
  z = factor(sample(letters[1:3], varN, replace = TRUE))
)
attrShorthandZ <- attr(
  dbarts::dbarts(
    y ~ x1 + z,
    data.frame(y = varY, dfFactorVar),
    control = varControl,
    variance = ~z,
    factors = "indicators"
  )$control,
  "bartcore.variance"
)
attrObjectZ <- attr(
  dbarts::dbarts(
    y ~ x1 + z,
    data.frame(y = varY, dfFactorVar),
    control = varControl,
    variance = dbarts::varianceForest(vars = ~z),
    factors = "indicators"
  )$control,
  "bartcore.variance"
)
expect_identical(attrShorthandZ, attrObjectZ)
expect_equal(attrShorthandZ$columns, 2:4) # z's 3 indicator columns

# variance = NULL/FALSE still declare no variance forest (the object's own
# NULL vars means "every column", the opposite of variance = NULL's meaning,
# by construction - only reachable through the object)
expect_null(varianceAttr(NULL))
expect_null(varianceAttr(FALSE))

# collision rule (3): the three removed formal spellings are now simply
# unknown arguments - bart2's dots rejection names each, unsuggested (none of
# the surviving formal names is close enough for agrep's default distance)
expect_error(
  fit2(y.gaussian, n.trees.variance = 10L),
  pattern = "^unknown argument 'n.trees.variance'$"
)
expect_error(
  fit2(y.gaussian, power.variance = 1),
  pattern = "^unknown argument 'power.variance'$"
)
expect_error(
  fit2(y.gaussian, base.variance = 1),
  pattern = "^unknown argument 'base.variance'$"
)
# dbarts()/dbartsSpec() have no dots at all, so R's own unused-argument
# refusal fires instead
expect_error(
  dbarts::dbarts(varX, varY, control = varControl, n.trees.variance = 10L),
  pattern = "unused argument"
)

# print/format smoke
printed <- capture.output(print(dbarts::varianceForest(n.trees = 20L)))
expect_true(any(grepl("variance forest", printed)))
expect_true(any(grepl("n.trees = 20", printed, fixed = TRUE)))
expect_true(any(grepl("<all columns>", printed, fixed = TRUE)))
formatted <- format(
  dbarts::varianceForest(vars = ~x1, n.trees = 20L, base = 0.9, power = 1.5)
)
expect_true(is.character(formatted))

# S7 (docs/plans/bart2-argument-consolidation.md 4.1, 4.2, 6.8): tree.prior/
# node.prior/resid.prior are now bart2 formals (NULL, appended after
# breaks/max.rows). A supplied object forwards unevaluated - exactly as k
# already does - and a shorthand that would otherwise build the same prior is
# a collision, refused by name; no object leaves the flat shorthand build
# untouched.

# reachability: node.prior = linear()/gp() and resid.prior = fixed() were
# unreachable from bart2 before S7 (no route existed to hand dbarts() a
# prior OBJECT). samplerOnly = TRUE returns bart2's sampler before any tree
# initialization, so it is a fair byte-identity comparison against a fresh
# dbarts() sampler built with the same object and control at the same seed -
# the real G3-closure evidence, not merely "it does not error".
reachSeed <- 909L
dfLeaf <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
yLeaf <- rnorm(n)

samplerViaBart2Linear <- dbarts::bart2(
  yLeaf ~ x1 + x2 + x3,
  dfLeaf,
  node.prior = dbarts::dbartsPriors$linear("x2"),
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = reachSeed,
  verbose = FALSE,
  samplerOnly = TRUE
)
samplerViaDbartsLinear <- dbarts::dbarts(
  yLeaf ~ x1 + x2 + x3,
  dfLeaf,
  node.prior = dbarts::dbartsPriors$linear("x2"),
  control = dbarts::dbartsControl(
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    seed = reachSeed,
    verbose = FALSE
  )
)
samplesViaBart2Linear <- samplerViaBart2Linear$run()
samplesViaDbartsLinear <- samplerViaDbartsLinear$run()
expect_identical(samplesViaBart2Linear$train, samplesViaDbartsLinear$train)
expect_identical(samplesViaBart2Linear$sigma, samplesViaDbartsLinear$sigma)

samplerViaBart2Fixed <- dbarts::bart2(
  x,
  y.gaussian,
  resid.prior = dbarts::dbartsPriors$fixed(1),
  n.trees = 3L,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = reachSeed,
  verbose = FALSE,
  samplerOnly = TRUE
)
samplerViaDbartsFixed <- dbarts::dbarts(
  x,
  y.gaussian,
  resid.prior = dbarts::dbartsPriors$fixed(1),
  control = dbarts::dbartsControl(
    n.trees = 3L,
    n.samples = 5L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    seed = reachSeed,
    verbose = FALSE
  )
)
samplesViaBart2Fixed <- samplerViaBart2Fixed$run()
samplesViaDbartsFixed <- samplerViaDbartsFixed$run()
expect_identical(samplesViaBart2Fixed$train, samplesViaDbartsFixed$train)
expect_identical(samplesViaBart2Fixed$sigma, samplesViaDbartsFixed$sigma)
# fixed(1) holds the residual VARIANCE at 1, not sigma bit-exactly at 1 (the
# rescale to/from the response's internal [-0.5, 0.5] coding is floating
# point), so the sanity check is a tolerance, not '=='
expect_true(all(abs(samplesViaBart2Fixed$sigma - 1) < 1e-8))

# collision refusals (4.2): a prior object plus a shorthand that would build
# it errors naming both, matching the dart/split.probs precedent's shape
expect_error(
  fit2(y.gaussian, tree.prior = dbarts::dbartsPriors$cgm(), power = 3),
  pattern = "'tree.prior' cannot be combined with 'power'"
)
expect_error(
  fit2(y.gaussian, node.prior = dbarts::dbartsPriors$normal(), k = 3),
  pattern = "'node.prior' cannot be combined with 'k'"
)
expect_error(
  fit2(y.gaussian, resid.prior = dbarts::dbartsPriors$fixed(1), sigdf = 5),
  pattern = "'resid.prior' cannot be combined with 'sigdf'"
)
# the m6 disposition: sigest collides too, even though it is not built into
# resid.prior by buildSamplerPriors (it rides dbarts()'s separate 'sigma')
expect_error(
  fit2(y.gaussian, resid.prior = dbarts::dbartsPriors$fixed(1), sigest = 2),
  pattern = "'resid.prior' cannot be combined with 'sigest'"
)

# family gating (3.c.2's new row): resid.prior joins the sigest/sigdf/
# sigquant trio - inert (silently overwritten with fixed(1), R/spec.R) under
# a fixed-unit-scale family, now diagnosed instead of silent
expect_warning(
  fit2(
    y.binary,
    family = "probit",
    resid.prior = dbarts::dbartsPriors$fixed(2)
  ),
  pattern = "resid.prior",
  class = "dbartsFamilyGatedWarning"
)
expect_equal(
  countWarnings(
    fit2(
      y.binary,
      family = "probit",
      resid.prior = dbarts::dbartsPriors$fixed(2)
    ),
    "dbartsFamilyGatedWarning"
  ),
  1L
)

# abbreviation breaks (4.1's table, 6.8): 'tree.prior' now collides with
# 'test' on 't=', and 'resid.prior' with 'resid.dist' on 'resid.='; R's own
# ambiguous partial-match error fires, not a package one
expect_error(
  dbarts::bart2(
    x,
    y.gaussian,
    t = 3L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "matches multiple formal arguments"
)
expect_error(
  dbarts::bart2(
    x,
    y.gaussian,
    resid. = 1,
    n.trees = 3L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "matches multiple formal arguments"
)

# tree.prior/node.prior are live on BOTH hurdle components, unlike the
# resid.prior/sigest quartet (positive half only): a differing tree.prior
# changes both components' draws relative to the default (a cheap liveness
# check, not a value-level assertion)
hurdleDefaultTree <- fit2(y.hurdle, family = "hurdle.lognormal", seed = 55L)
hurdleOtherTree <- fit2(
  y.hurdle,
  family = "hurdle.lognormal",
  tree.prior = dbarts::dbartsPriors$cgm(power = 4),
  seed = 55L
)
expect_false(identical(
  hurdleDefaultTree$occupancy$yhat.train,
  hurdleOtherTree$occupancy$yhat.train
))
expect_false(identical(
  hurdleDefaultTree$positive$yhat.train,
  hurdleOtherTree$positive$yhat.train
))
expect_true(any(grepl("base\\s*= 0.9", formatted)))

# S8 (docs/plans/bart2-argument-consolidation.md 3.f, 4.3, 6.2): xbart loses
# control = entirely; n.cuts/useQuantiles/n.thin/storage/tree.prior are
# appended flat formals instead. xbart has no dots (fork 5's rejection-only
# channel never reached it), so control = is a native unused-argument error.
xbartKnobs <- c("n.cuts", "useQuantiles", "n.thin", "storage")
expect_true(setequal(
  intersect(xbartKnobs, names(formals(dbarts::xbart))),
  xbartKnobs
))
for (knob in xbartKnobs) {
  expect_equal(
    deparse(formals(dbarts::xbart)[[knob]]),
    deparse(formals(dbarts::dbartsControl)[[knob]])
  )
}
expect_false("control" %in% names(formals(dbarts::xbart)))
expect_false("..." %in% names(formals(dbarts::xbart)))
expect_equal(length(formals(dbarts::xbart)), 32L)
expect_error(
  dbarts::xbart(
    x,
    y.gaussian,
    control = dbarts::dbartsControl(),
    n.reps = 1L,
    n.threads = 1L
  ),
  pattern = "unused argument"
)
