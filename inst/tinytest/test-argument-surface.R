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
# tree.prior/node.prior/resid.prior are not yet bart2 formals (S7); their
# rows are inert until then, since the loop only walks names shared today.
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
