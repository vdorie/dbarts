# Public Bayesian causal forest creation (docs/design/bcf.md): dbarts() and
# dbartsSpec() accept a 0/1 treatment vector and build an ordinary
# dbartsSampler holding the two-forest model. The internal
# dbarts:::bartcoreBCFSampler route is the oracle - identical (control, model,
# data) and a fixed rng seed must reproduce it draw for draw - and every option
# the two-forest chain does not read must refuse at creation rather than be
# dropped in silence.

set.seed(3)
n <- 300L
p <- 4L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

seededControl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 50L,
    n.samples = 10L,
    updateState = FALSE,
    rngSeed = 17L,
    ...
  )
}

# a public sampler's raw handle, so the internal per-forest readers apply to it
handleOf <- function(sampler) list(ptr = sampler$getPointer())

# --- F1, positive half: with control@rngSeed set every chain's generator is
# independent of R's stream, so the public path and the internal one receive
# identical seeds from identical specifications and must agree bitwise on all
# six channels the bcf-equivalence fixture reports. ---
control <- seededControl()
publicSampler <- dbarts(
  x,
  y,
  treatment = z,
  treatmentForest = treatmentForest(n.trees = 25L, sd.moderate = 1.5),
  control = control
)
publicResult <- publicSampler$run(0L, 10L)

internalHost <- dbarts(x, y, control = control)
internalSampler <- dbarts:::bartcoreBCFSampler(
  internalHost,
  z,
  n.trees.treatment = 25L,
  sd.moderate = 1.5
)
internalResult <- dbarts:::bartcoreRun(internalSampler, 0L, 10L)

expect_identical(publicResult$train, internalResult$train)
expect_identical(publicResult$sigma, internalResult$sigma)
expect_identical(publicResult$varcount, internalResult$varcount)
expect_identical(
  dbarts:::bartcoreForestFits(handleOf(publicSampler), 0L),
  dbarts:::bartcoreForestFits(internalSampler, 0L)
)
expect_identical(
  dbarts:::bartcoreForestFits(handleOf(publicSampler), 1L),
  dbarts:::bartcoreForestFits(internalSampler, 1L)
)
expect_identical(
  dbarts:::bartcoreBCFGlue(handleOf(publicSampler)),
  dbarts:::bartcoreBCFGlue(internalSampler)
)

# the sampler this produced is an ordinary dbartsSampler carrying two forests
expect_true(inherits(publicSampler, "dbartsSampler"))
expect_identical(publicSampler$data@treatment, as.double(z))
expect_true(all(is.finite(publicResult$train)))
expect_true(all(publicResult$sigma > 0))

# --- F1, negative half: UNSEEDED, the internal route builds and discards a
# host engine first, which draws n.chains unif_rand()s off R's stream, so the
# two routes must NOT agree. Written explicitly so a divergence here is read as
# the expected stream offset rather than as a creation bug. ---
unseeded <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 50L,
  n.samples = 10L,
  updateState = FALSE
)
set.seed(99L)
unseededPublic <- dbarts(x, y, treatment = z, control = unseeded)$run(0L, 5L)
set.seed(99L)
unseededInternal <- dbarts:::bartcoreRun(
  dbarts:::bartcoreBCFSampler(dbarts(x, y, control = unseeded), z),
  0L,
  5L
)
expect_false(identical(unseededPublic$train, unseededInternal$train))

# --- the reported draws really are the two-forest blend: the glue and the
# per-forest fits reconstruct the recorded train draw through the stored
# response transform (the S0 pin's identity, now on a public sampler) ---
glue <- dbarts:::bartcoreBCFGlue(handleOf(publicSampler))
muFits <- dbarts:::bartcoreForestFits(handleOf(publicSampler), 0L)[, 1L]
tauFits <- dbarts:::bartcoreForestFits(handleOf(publicSampler), 1L)[, 1L]
fitScale <- dbarts:::bartcoreStoreState(handleOf(publicSampler))[[1L]]$fit.scale
scale <- fitScale[2L] - fitScale[1L]
shift <- scale * 0.5 + fitScale[1L]
bz <- ifelse(z != 0, glue[3L, 1L], glue[2L, 1L])
expect_equal(
  scale * (glue[1L, 1L] * muFits + bz * tauFits) + shift,
  publicResult$train[, 10L, 1L],
  tolerance = 1e-10
)

# --- the moderator restriction reaches the treatment forest: tau splits only
# on the declared columns, mu on whatever it likes ---
restricted <- dbarts(
  x,
  y,
  treatment = z,
  moderators = c("x3", "x4"),
  control = seededControl()
)
restricted$run(0L, 10L)
tauCounts <- dbarts:::bartcoreForestVariableCounts(handleOf(restricted), 1L)
expect_true(all(tauCounts[1:2, ] == 0L))
expect_true(sum(tauCounts[3:4, ]) > 0L)

# --- dbartsSpec() reaches the same model, and carries the treatment on the
# data object rather than on the control ---
specSampler <- do.call(
  function(control, model, data, family) {
    new("dbartsSampler", control, model, data)
  },
  dbartsSpec(
    dbartsData(x, y),
    seededControl(),
    treatment = z,
    treatmentForest = treatmentForest(n.trees = 25L, sd.moderate = 1.5)
  )
)
expect_identical(specSampler$run(0L, 10L)$train, publicResult$train)

# a treatment supplied to dbartsData() rides the data object and is restricted
# by 'subset' exactly as the weights beside it are
subsetData <- dbartsData(x, y, subset = 1:100, treatment = z)
expect_identical(subsetData@treatment, as.double(z[1:100]))
expect_error(dbartsData(x, y, treatment = rep(2, n)), "coded 0")
expect_error(dbartsData(x, y, treatment = z[1:10]), "length of 'treatment'")

# --- F4: every option the two-forest chain does not read refuses at creation,
# one assertion each. A silently accepted one is the failure mode this
# creation route exists to prevent. ---
expect_error(
  dbarts(x, y, treatment = z, tree.prior = dart, control = seededControl()),
  "DART tree prior"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    tree.prior = cgm(split.probs = c(0.4, 0.2, 0.2, 0.2)),
    control = seededControl()
  ),
  "split.probs"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    monotone = c(1, 0, 0, 0),
    control = seededControl()
  ),
  "monotone"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    node.prior = linear(columns = "x1"),
    control = seededControl()
  ),
  "linear node prior"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    node.prior = gp(columns = "x1"),
    control = seededControl()
  ),
  "Gaussian-process node prior"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    node.prior = normal(chi(1.5)),
    control = seededControl()
  ),
  "'k' hyperprior"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    node.prior = normal(3),
    control = seededControl()
  ),
  "non-default 'k'"
)
expect_error(
  dbarts(x, y, treatment = z, variance = TRUE, control = seededControl()),
  "'variance'"
)
expect_error(
  dbarts(x, y, treatment = z, control = seededControl(storage = "single")),
  "storage = \"single\""
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0.1,
      change = 0.3,
      birth = 0.5
    ),
    control = seededControl()
  ),
  "non-default 'proposal.probs'"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    resid.dist = student(5),
    control = seededControl()
  ),
  "Student-t residuals"
)
expect_error(
  dbarts(x, y, test = x, treatment = z, control = seededControl()),
  "test predictors"
)
# per-column cut counts survive only through dbartsSpec(), which keeps a data
# object's own resolved n.cuts
unevenCuts <- dbartsData(x, y, treatment = z)
unevenCuts@n.cuts <- c(50L, 100L, 100L, 100L)
expect_error(dbartsSpec(unevenCuts, seededControl()), "per-column 'n.cuts'")
# grouped random effects reach creation on rbart_vi's internal control
# attribute; the composition with a second forest is a door, not a build
groupedControl <- seededControl()
attr(groupedControl, "bartcore.groups") <- list(
  indices = rep(1L, n),
  n.groups = 1L,
  prior = "cauchy",
  rel.scale = 1,
  n.steps = 1L
)
expect_error(
  dbartsSpec(dbartsData(x, y, treatment = z), groupedControl),
  "grouped random effects"
)
# a latent-scale family owns the precision channel the glue routes through
expect_error(
  dbarts(x, rbinom(n, 1L, 0.5), treatment = z, control = seededControl()),
  "gaussian"
)
# the two companion arguments configure a forest a treatment vector selects
expect_error(
  dbarts(x, y, moderators = "x3", control = seededControl()),
  "'treatment' vector selects"
)
expect_error(
  dbarts(x, y, treatmentForest = treatmentForest(), control = seededControl()),
  "'treatment' vector selects"
)

# --- F5: a creation the factory would refuse raises an R condition, and no
# usable handle escapes. Driven past the R-layer refusal by hand-attaching the
# variance attribute to an already-resolved BCF spec, so the bridge's own
# backstop is what answers. ---
resolved <- dbartsSpec(dbartsData(x, y), seededControl(), treatment = z)
withVariance <- resolved
attr(withVariance$control, "bartcore.variance") <- list(
  n.trees = 20L,
  base = 0.95,
  power = 2.0,
  columns = NULL
)
escapedSampler <- NULL
expect_error(
  escapedSampler <- new(
    "dbartsSampler",
    withVariance$control,
    withVariance$model,
    withVariance$data
  ),
  "variance forest"
)
# the condition is raised before any handle is bound, so no external pointer -
# least of all one wrapping a null sampler - escapes to be used
expect_null(escapedSampler)

# the bridge backstop also answers a hand-built model the R layer never emits
scaledModel <- resolved
scaledModel$model@node.scale <- 3.0
expect_error(
  new(
    "dbartsSampler",
    scaledModel$control,
    scaledModel$model,
    scaledModel$data
  ),
  "non-default node scale"
)

# --- F7: the two halves of the specification are cross-checked in BOTH
# directions, so a stripped one is a loud refusal naming the missing piece
# rather than a silent single-forest fit ---
strippedConfig <- resolved
attr(strippedConfig$control, "bartcore.bcf") <- NULL
expect_error(
  new(
    "dbartsSampler",
    strippedConfig$control,
    strippedConfig$model,
    strippedConfig$data
  ),
  "no treatment forest was configured"
)
strippedTreatment <- resolved
strippedTreatment$data@treatment <- NULL
expect_error(
  new(
    "dbartsSampler",
    strippedTreatment$control,
    strippedTreatment$model,
    strippedTreatment$data
  ),
  "carry no treatment vector"
)

# --- treatmentForest() validates its own knobs at fit time ---
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    treatmentForest = treatmentForest(base = 1.5),
    control = seededControl()
  ),
  "'base'"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    treatmentForest = treatmentForest(n.trees = 0L),
    control = seededControl()
  ),
  "'n.trees'"
)
expect_error(
  dbarts(
    x,
    y,
    treatment = z,
    treatmentForest = list(n.trees = 25L),
    control = seededControl()
  ),
  "treatmentForest\\(\\) specification"
)
