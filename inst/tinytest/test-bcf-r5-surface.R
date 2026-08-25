# The dbartsSampler surface over a public Bayesian causal forest:
# $setForestBasis mirrors the engine
# and data@bases, $getForestFits/$getForestAmplitudes/
# $getForestVariableCounts read the per-forest channels the low-level
# bartcoreForestFits/bartcoreForestAmplitudes
# route already exposed, refused mutations get a BCF-specific message rather
# than the generic multi-forest one, and $setControl carries bartcore.*
# control attributes forward so a save/load round trip re-creates cleanly.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(5)
n <- 200L
p <- 4L
x <- matrix(runif(n * p), n, p)
z <- rbinom(n, 1L, 0.5)
mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
tau <- 1 + 2 * x[, 3L]
y <- mu + z * tau + rnorm(n, sd = 0.2)

seededControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 40L,
    n.samples = 10L,
    updateState = FALSE,
    seed = 23L,
    ...
  )
}

# --- getForestFits/getForestAmplitudes/getForestVariableCounts read the same
# the low-level route does, and the driver-loop reconstruction identity
# holds through the R-level accessors too ---
sampler <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
result <- sampler$run(0L, 5L)

muFits <- sampler$getForestFits(1L)
tauFits <- sampler$getForestFits(2L)
glue <- sampler$getForestAmplitudes()
muCounts <- sampler$getForestVariableCounts(1L)
tauCounts <- sampler$getForestVariableCounts(2L)

expect_equal(dim(muFits), c(n, 1L))
expect_equal(dim(tauFits), c(n, 1L))
expect_equal(dim(glue), c(3L, 1L))
expect_equal(dim(muCounts), c(p, 1L))
expect_equal(dim(tauCounts), c(p, 1L))
expect_true(all(is.finite(muFits)) && all(is.finite(tauFits)))

# identical to the same low-level readers on the same live pointer: the
# R-level methods add no computation of their own, and their forest index is 1-based
# (1 = prognostic, 2 = treatment) against the low-level route's 0-based one
lowLevel <- list(ptr = sampler$getPointer())
expect_identical(muFits, bartcoreForestFits(lowLevel, 0L))
expect_identical(tauFits, bartcoreForestFits(lowLevel, 1L))
expect_identical(glue, bartcoreForestAmplitudes(lowLevel))
expect_identical(
  muCounts,
  bartcoreForestVariableCounts(lowLevel, 0L)
)
expect_identical(
  tauCounts,
  bartcoreForestVariableCounts(lowLevel, 1L)
)

# 0 is refused rather than silently naming the prognostic forest
# (resolveForestIndex, shared with setForestWeights/getCalibration)
expect_error(sampler$getForestFits(0L), "single positive integer")
expect_error(sampler$getForestVariableCounts(0L), "single positive integer")

# getForestVariableCounts names its rows from data@x's colnames when present
# (x above carries none, the no-op guard) and leaves the values themselves
# unchanged from the unnamed path
expect_null(dimnames(muCounts))
namedX <- x
colnames(namedX) <- paste0("v", seq_len(p))
namedSampler <- dbarts(
  namedX,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
namedSampler$run(0L, 5L)
namedCounts <- namedSampler$getForestVariableCounts(1L)
expect_identical(rownames(namedCounts), colnames(namedX))
expect_identical(unname(namedCounts), muCounts)

sampler$storeState()
fitScale <- sampler$state[[1L]]$fit.scale
scale <- fitScale[2L] - fitScale[1L]
shift <- scale * 0.5 + fitScale[1L]
bz <- ifelse(z != 0, glue[3L, 1L], glue[2L, 1L])
expect_equal(
  scale * (glue[1L, 1L] * muFits[, 1L] + bz * tauFits[, 1L]) + shift,
  result$train[, 5L],
  tolerance = 1e-10
)

# --- BCF-specific messages on the refused mutations: an R-level BCF sampler names
# the amplitude capability rather than surfacing the bridge's generic
# "multi-forest" wording ---
refused <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)

expect_error(refused$setData(refused$data), "carries forest amplitudes")
expect_error(refused$setModel(refused$model), "carries forest amplitudes")
expect_error(
  refused$setResponse(y, updateScale = TRUE),
  "carries forest amplitudes"
)
expect_error(
  refused$setOffset(rep(0, n), updateScale = TRUE),
  "carries forest amplitudes"
)
# the transactional column update is no longer refused: revalidateAllChains
# loops both forests, so it installs under the empty-leaf veto and reports a
# logical rather than erroring. Replacing a column with its own values cannot
# empty a leaf, so the answer is TRUE
expect_true(refused$setPredictor(x[, 1L], column = 1L))
# the per-observation session is no longer refused either: its cell guard
# caches every forest, pruned to the trees the column can move, so it returns
# the install mask rather than erroring. Re-installing the column's own values
# moves no observation, so every row installs
installed.r5 <- refused$setPredictor(
  x[, 1L],
  column = 1L,
  forceUpdate = "partial"
)
expect_true(is.logical(installed.r5) && all(installed.r5))
# a two-level replacement returns a mask of the same shape and leaves the
# sampler runnable. This sampler has not been run, so its trees are stumps and
# no row can empty a leaf; the DECLINE half of the veto is pinned on the
# burned-in low-level handles (test-bcf-mutation-pins.R,
# test-multi-forest-seam.R), where a leaf can be down to its last occupant
installed.coarse <- refused$setPredictor(
  ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
  column = 1L,
  forceUpdate = "partial"
)
expect_true(is.logical(installed.coarse) && length(installed.coarse) == n)
expect_true(all(is.finite(refused$run(0L, 5L)$train)))

# the accepted BCF mutations are unaffected by the new pre-checks
expect_silent(refused$setPredictor(x, forceUpdate = TRUE))
expect_silent(refused$setResponse(y, updateScale = FALSE))
expect_silent(refused$setOffset(rep(0, n), updateScale = FALSE))
expect_silent(refused$setForestBasis(2L, factor(z)))
# any forest takes a basis, the first included - it widens that forest's
# amplitude block from the implicit intercept's single coordinate to one per
# level - while the index itself still refuses on both sides
expect_silent(refused$setForestBasis(1L, factor(z)))
expect_equal(nrow(refused$getForestAmplitudes(1L)), 2L)
expect_error(refused$setForestBasis(0L, factor(z)), "single positive integer")
expect_error(refused$setForestBasis(3L, factor(z)), "out of range")

# the test surface's own refusal (refuseUndefinedTestFits); light regression
# coverage on the R-level path, previously pinned only on the low-level
# handle. setTestOffset is not exercised here: with no test matrix (this
# sampler can never install one - setTestPredictor is refused above), the
# R-level method's own "test matrix is NULL" precondition fires before the
# .Call, so the bridge's message is unreachable through this method - a
# pre-existing, capability-independent check untouched here.
expect_error(
  refused$predict(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)
expect_error(
  refused$setTestPredictor(x[1:5, , drop = FALSE]),
  "have no off-sample basis"
)

# regression canary: the new pre-checks fire on the capability alone and must
# not reach an
# ordinary single-forest sampler (exercised at scale by the rest of the suite)
plain <- dbarts(x, y, control = seededControl())
expect_silent(plain$setResponse(y, updateScale = TRUE))

# --- the basis-mirror contract, both halves: $setForestBasis mirrors the
# engine AND data@bases, and a save/load round trip continues from the
# mirrored assignment rather than the one the sampler was created with ---
mirror <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
mirror$run(10L, 5L)
z2 <- rep(0, n)
# an explicit level set, so the basis still has the two columns the amplitudes
# ride on when the data happen to sit entirely in one of them
mirror$setForestBasis(2L, factor(z2, levels = c(0, 1)))
expect_equal(mirror$data@bases[[2L]], unname(cbind(1 - z2, z2)))

mirror$storeState()
mirrorFile <- tempfile(fileext = ".rds")
saveRDS(mirror, mirrorFile)
reloadedMirror <- readRDS(mirrorFile)
unlink(mirrorFile)

# forces getPointer()'s transparent re-creation (the external pointer does
# not survive serialization); the slot rides the object, so this needs no
# extra plumbing
reloadedResult <- reloadedMirror$run(0L, 1L)
reloadedGlue <- reloadedMirror$getForestAmplitudes()
reloadedMuFits <- reloadedMirror$getForestFits(1L)[, 1L]
reloadedTauFits <- reloadedMirror$getForestFits(2L)[, 1L]

# scale/shift (fit.scale = range(y)) are reused from the sampler above: y is
# identical across every sampler in this file, so the affine map back to the
# reported scale does not depend on which one built it
reconstructTrain <- function(zVec) {
  bzVec <- ifelse(zVec != 0, reloadedGlue[3L, 1L], reloadedGlue[2L, 1L])
  scale *
    (reloadedGlue[1L, 1L] * reloadedMuFits + bzVec * reloadedTauFits) +
    shift
}

# POSITIVE: reconstructing with z2 (the mirrored assignment) matches
diffZ2 <- max(abs(reconstructTrain(z2) - reloadedResult$train[, 1L]))
expect_true(diffZ2 < 1e-4)
# NEGATIVE: reconstructing with the creation-time z does not - the mirror is
# what is proven here, not merely that re-creation runs at all. b0 and b1 are
# continuous Gibbs draws and tau(x) is a real fitted forest, so this margin
# only requires b0 != b1 to double precision - not a particular separation -
# which holds with probability 1 for any non-degenerate run
diffOriginalZ <- max(abs(reconstructTrain(z) - reloadedResult$train[, 1L]))
expect_true(diffOriginalZ > 1e-6)

# --- BCF leg: setControl preserves attr(control, "bartcore.forests"), so
# getPointer()'s re-creation after a save/load round trip succeeds and
# carries the same bases. Pre-fix, this raised "the data carry forest bases
# but no basis forest was configured": setControl replaced control wholesale,
# dropping the attribute, so createHolder saw data@bases with no matching
# forest configuration (the bidirectional cross-check firing on the
# half setControl silently orphaned). ---
controlled <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = seededControl()
)
controlled$run(0L, 5L)
controlled$storeState()
controlled$setControl(seededControl(printEvery = 50L))

controlledFile <- tempfile(fileext = ".rds")
saveRDS(controlled, controlledFile)
reloadedControlled <- readRDS(controlledFile)
unlink(controlledFile)

reloadedControlled$getPointer()
expect_equal(
  reloadedControlled$data@bases[[2L]],
  cbind(1 - as.double(z), as.double(z))
)
rerunControlled <- reloadedControlled$run(0L, 1L)
expect_true(all(is.finite(rerunControlled$train)))

# --- heteroscedastic leg: the same fix, on the model type whose defect
# shows up as "state is not consistent with this sampler" - setState's
# forest-block-count check, since a variance forest has no data-side vector
# to cross-check against, so the mismatch was silent until the state, not
# the creation, was pushed ---
varControl <- seededControl()
varSampler <- dbarts(
  x,
  y,
  variance = varianceForest(n.trees = 10L),
  control = varControl
)
varSampler$run(0L, 5L)
varSampler$storeState()
varSampler$setControl(seededControl(printEvery = 50L))

varFile <- tempfile(fileext = ".rds")
saveRDS(varSampler, varFile)
reloadedVar <- readRDS(varFile)
unlink(varFile)

reloadedVar$getPointer()
reloadedVarResult <- reloadedVar$run(0L, 1L)
expect_true(!is.null(reloadedVarResult$variance))
expect_equal(length(reloadedVarResult$variance), n)
expect_true(all(is.finite(reloadedVarResult$variance)))
expect_true(all(reloadedVarResult$variance > 0))
