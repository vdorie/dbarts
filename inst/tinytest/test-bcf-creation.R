# Public multi-forest creation (docs/design/bcf.md): dbarts() and dbartsSpec()
# take a forests = list(forest(...)) declaration and build an ordinary
# dbartsSampler holding the model it names. Two forests, the second carrying a
# two-level factor basis, is the Bayesian causal forest. The internal
# dbarts:::bartcoreBCFSampler route is the oracle - identical (control, model,
# data) and a fixed rng seed must reproduce it draw for draw - and every option
# the two-forest chain does not read, and every declaration today's engine
# cannot honour, must refuse at creation rather than be dropped in silence.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

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
    seed = 17L,
    ...
  )
}

# a public sampler's raw handle, so the internal per-forest readers apply to it
handleOf <- function(sampler) list(ptr = sampler$getPointer())

# the declaration every refusal below is attached to: a plain first forest and
# a second one whose two-level factor basis carries the (b0, b1) amplitudes
twoForests <- list(forest(), forest(basis = ~ factor(z)))

# --- the creation-reproduction contract, positive half: with control@seed
# set every chain's generator is
# independent of R's stream, so the public path and the internal one receive
# identical seeds from identical specifications and must agree bitwise on all
# six channels the bcf-equivalence fixture reports. The removed treatment =
# route needs no reconstruction: the internal constructor is, and always was,
# the oracle this pin compares against. ---
control <- seededControl()
publicSampler <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z), n.trees = 25L, sd = 1.5)
  ),
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
internalResult <- bartcoreRun(internalSampler, 0L, 10L)

expect_identical(publicResult$train, internalResult$train)
expect_identical(publicResult$sigma, internalResult$sigma)
expect_identical(publicResult$varcount, internalResult$varcount)
expect_identical(
  bartcoreForestFits(handleOf(publicSampler), 0L),
  bartcoreForestFits(internalSampler, 0L)
)
expect_identical(
  bartcoreForestFits(handleOf(publicSampler), 1L),
  bartcoreForestFits(internalSampler, 1L)
)
expect_identical(
  bartcoreForestAmplitudes(handleOf(publicSampler)),
  bartcoreForestAmplitudes(internalSampler)
)

# the sampler this produced is an ordinary dbartsSampler carrying two forests,
# and the level indicators the factor basis expanded to are the ones the slot
# takes, forest by forest - the first forest declaring none
expect_true(inherits(publicSampler, "dbartsSampler"))
expect_null(publicSampler$data@bases[[1L]])
expect_equal(
  publicSampler$data@bases[[2L]],
  cbind(1 - as.double(z), as.double(z))
)
expect_true(all(is.finite(publicResult$train)))
expect_true(all(publicResult$sigma > 0))

# --- the creation-reproduction contract, negative half, arm 1: the
# expansion's level ORDER is load-bearing.
# Swapping the two indicator columns swaps which rows b1 scales, so a sampler
# built from the reversed factor must NOT reproduce the oracle. Without this a
# silently transposed basis would pass the positive half on symmetry alone. ---
swappedSampler <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z, levels = c(1, 0)), n.trees = 25L, sd = 1.5)
  ),
  control = seededControl()
)
expect_equal(
  swappedSampler$data@bases[[2L]],
  cbind(as.double(z), 1 - as.double(z))
)
expect_false(identical(
  swappedSampler$run(0L, 10L)$train,
  internalResult$train
))

# --- the creation-reproduction contract, negative half, arm 2: UNSEEDED,
# the internal route builds and
# discards a host engine first, which draws n.chains unif_rand()s off R's
# stream, so the two routes must NOT agree. Written explicitly so a divergence
# here is read as the expected stream offset rather than as a creation bug. ---
unseeded <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 50L,
  n.samples = 10L,
  updateState = FALSE
)
set.seed(99L)
unseededPublic <- dbarts(x, y, forests = twoForests, control = unseeded)$run(
  0L,
  5L
)
set.seed(99L)
unseededInternal <- bartcoreRun(
  dbarts:::bartcoreBCFSampler(dbarts(x, y, control = unseeded), z),
  0L,
  5L
)
expect_false(identical(unseededPublic$train, unseededInternal$train))

# --- the reported draws really are the two-forest blend: the amplitudes and
# the per-forest fits reconstruct the recorded train draw through the stored
# response transform (the identity a low-level pin already established, now
# on a public sampler) ---
glue <- bartcoreForestAmplitudes(handleOf(publicSampler))
muFits <- bartcoreForestFits(handleOf(publicSampler), 0L)[, 1L]
tauFits <- bartcoreForestFits(handleOf(publicSampler), 1L)[, 1L]
fitScale <- bartcoreStoreState(handleOf(publicSampler))[[1L]]$fit.scale
scale <- fitScale[2L] - fitScale[1L]
shift <- scale * 0.5 + fitScale[1L]
bz <- ifelse(z != 0, glue[3L, 1L], glue[2L, 1L])
expect_equal(
  scale * (glue[1L, 1L] * muFits + bz * tauFits) + shift,
  publicResult$train[, 10L, 1L],
  tolerance = 1e-10
)

# --- the second forest's vars restriction reaches it: tau splits only on the
# declared columns, mu on whatever it likes ---
restricted <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), vars = c("x3", "x4"))),
  control = seededControl()
)
restricted$run(0L, 10L)
tauCounts <- bartcoreForestVariableCounts(handleOf(restricted), 1L)
expect_true(all(tauCounts[1:2, ] == 0L))
expect_true(sum(tauCounts[3:4, ]) > 0L)

# --- the whole knob map is honoured, not accepted and dropped. Each knob of
# the second forest reaches the eight doubles the control attribute carries,
# and the FIRST forest's structural knobs restate the fit's own control and
# tree prior rather than adding a second set. A knob that went nowhere would
# leave a default here. ---
knobs <- dbartsSpec(
  dbartsData(x, y),
  seededControl(),
  forests = list(
    forest(
      n.trees = 30L,
      base = 0.7,
      power = 1.5,
      sd = 2.5,
      update.amplitude = FALSE
    ),
    forest(
      basis = ~ factor(z),
      n.trees = 15L,
      base = 0.4,
      power = 2,
      sd = 1.25,
      amplitude.prior.variance = 0.75,
      update.amplitude = FALSE
    )
  )
)
# ragged: one length-8 vector per forest. The first forest declares no basis,
# so its `sd` is the half-Cauchy median of a plain scalar amplitude and its
# node scale stays at the response sd; the second declares one, so its `sd`
# rides the node scale through the half-normal median and its
# amplitude.prior.variance is the fixed prior on the block
expect_equal(
  attr(knobs$control, "bartcore.forests")$params,
  list(
    c(30, 0.7, 1.5, 1, 1, 1, 2.5, 0),
    c(15, 0.4, 2, 1.25, 0.674, 0.75, 0, 0)
  )
)
expect_identical(knobs$control@n.trees, 30L)
expect_equal(knobs$model@tree.prior@base, 0.7)
expect_equal(knobs$model@tree.prior@power, 1.5)

# --- a forest() n.trees/base/power on the first forest is a more specific
# declaration than an explicitly supplied control@n.trees or tree.prior of the
# same name, and governs the fit over it rather than being silently overridden
# by it: the run really is the seven-tree, base = 0.3, power = 1.1 fit, not
# the hundred-tree, base = 0.6, power = 4 one it disagrees with ---
disagreeingControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 100L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L
)
sevenTreeControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 7L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 7L, base = 0.3, power = 1.1)),
    tree.prior = cgm(base = 0.6, power = 4),
    control = disagreeingControl
  )$run(0L, 5L)$train,
  dbarts(
    x,
    y,
    tree.prior = cgm(base = 0.3, power = 1.1),
    control = sevenTreeControl
  )$run(0L, 5L)$train
)

# The pinned-amplitude public fit: update.amplitude = FALSE on both
# forests is the pinned-amplitude model (y = mu + z tau exactly) the on-ramp
# vignette's continuity falsifier composes against, and it must reproduce the
# internal route's update.a = update.b = FALSE draw for draw
pinnedPublic <- dbarts(
  x,
  y,
  forests = list(
    forest(update.amplitude = FALSE),
    forest(basis = ~ factor(z), n.trees = 25L, update.amplitude = FALSE)
  ),
  control = seededControl()
)
pinnedInternal <- dbarts:::bartcoreBCFSampler(
  dbarts(x, y, control = seededControl()),
  z,
  n.trees.treatment = 25L,
  update.a = FALSE,
  update.b = FALSE
)
expect_identical(
  pinnedPublic$run(0L, 10L)$train,
  bartcoreRun(pinnedInternal, 0L, 10L)$train
)

# --- forests = NULL is byte-neutral, and a single-forest declaration is
# the same fit with its structural knobs restated ---
plainResult <- dbarts(x, y, control = seededControl())$run(0L, 5L)
expect_identical(
  dbarts(x, y, forests = NULL, control = seededControl())$run(0L, 5L)$train,
  plainResult$train
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 50L)),
    control = seededControl()
  )$run(0L, 5L)$train,
  plainResult$train
)
# and a restated knob really is the fit's own: a different tree count draws
# what the control carrying that count draws
twentyTrees <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 20L)),
    control = seededControl()
  )$run(0L, 5L)$train,
  dbarts(x, y, control = twentyTrees)$run(0L, 5L)$train
)

# --- dbartsSpec() reaches the same model, and carries the basis column on the
# data object rather than on the control ---
specSampler <- do.call(
  function(control, model, data, family) {
    new("dbartsSampler", control, model, data)
  },
  dbartsSpec(
    dbartsData(x, y),
    seededControl(),
    forests = list(
      forest(),
      forest(basis = ~ factor(z), n.trees = 25L, sd = 1.5)
    )
  )
)
expect_identical(specSampler$run(0L, 10L)$train, publicResult$train)

# bases supplied to dbartsData() ride the data object and are restricted by
# 'subset' exactly as the weights beside them are; that argument is the
# treatment slot's successor, and it still selects the multi-forest model with
# no declaration at all
zBasis <- cbind(1 - as.double(z), as.double(z))
subsetData <- dbartsData(x, y, subset = 1:100, bases = list(NULL, zBasis))
expect_equal(subsetData@bases[[2L]], zBasis[1:100, ])
expect_error(
  dbartsData(x, y, bases = list(NULL, cbind(NA_real_, 1))),
  "must have the same length as 'y'"
)
expect_error(
  dbartsData(x, y, bases = list(NULL, zBasis[1:10, ])),
  "'bases' must have the same length"
)
# a non-finite value is refused rather than reaching the multiplier
expect_error(
  dbartsData(x, y, bases = list(NULL, cbind(rep(Inf, n), 0))),
  "must all be finite"
)

# ... while a caller who declared a forest's BASIS is refused in that word, on
# every surface that takes one, rather than in the data object's
expect_error(
  dbartsSpec(
    dbartsData(x, y),
    seededControl(),
    forests = list(forest(), forest(basis = ~ factor(z[1:10])))
  ),
  "'basis' must have the same length"
)
expect_error(
  dbarts(x, y, forests = list(forest(), forest(basis = ~ factor(z[1:10])))),
  "'basis' must have the same length"
)
expect_error(
  dbarts(
    x,
    y,
    control = seededControl(),
    forests = list(forest(), forest(basis = ~ factor(z)))
  )$setForestBasis(2L, ~ factor(z[1:10])),
  "'basis' must have the same length"
)
offSlot <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl()
)
# the transport is RAGGED - one length-8 vector per forest, not one length-8
# vector for the pair - and which magnitude channel `sd` reaches is decided by
# whether the forest carries a basis
expect_equal(
  attr(offSlot$control, "bartcore.forests")$params,
  list(c(50, 0.25, 3, 1, 1, 1, 2, 1), c(50, 0.25, 3, 1, 0.674, 0.5, 0, 1))
)
# the same transport under a LATENT family, which no assertion above has
# run: `sd` still reaches slot 4 on a forest carrying a basis and slot 7 on one
# carrying none, and the family reaches the model without disturbing either.
# The two declared values DIFFER, and neither is 1, so a confusion between the
# slots - or between a forest's own vector and its neighbour's - is visible;
# both defaults are 1 in this stretch of the vector and would not be.
latentSlots <- dbartsSpec(
  dbartsData(x, as.double(y > median(y)), bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(forest(sd = 3.5), forest(sd = 1.25)),
  family = "logistic"
)
expect_equal(latentSlots$model@family, "logistic")
expect_equal(
  attr(latentSlots$control, "bartcore.forests")$params,
  list(c(50, 0.25, 3, 1, 1, 1, 3.5, 1), c(50, 0.25, 3, 1.25, 0.674, 0.5, 0, 1))
)
# and a basis declared alongside subset reaches dbartsData's own alignment
subsetForests <- dbarts(
  x,
  y,
  subset = 1:100,
  forests = twoForests,
  control = seededControl()
)
expect_equal(subsetForests$data@bases[[2L]], zBasis[1:100, ])

# --- the hasBasis escape (model.R's validateForestKnobs/resolveForests): a
# basis reaching the data through dbartsData(bases=) rather
# than a forest() basis= unlocks the amplitude knobs and a basis-less second
# forest just as a declared basis would, and those knobs really reach the
# resolved spec ---
sdViaHasBasis <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(forest(sd = 3.5))
)
expect_equal(
  attr(sdViaHasBasis$control, "bartcore.forests")$params[[1L]][7L],
  3.5
)

noUpdateViaHasBasis <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(forest(update.amplitude = FALSE))
)
expect_equal(
  attr(noUpdateViaHasBasis$control, "bartcore.forests")$params[[1L]][8L],
  0
)

# a second forest with no basis of its own still takes its own knobs when the
# basis reaches the data another way
noBasisSecondForest <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(forest(), forest(n.trees = 25L))
)
expect_equal(
  attr(noBasisSecondForest$control, "bartcore.forests")$params[[2L]][1L],
  25
)

# --- interactions()/blocks() declared on a forest() reach the resolved
# per-forest constraint the control attribute carries, not only the top-level
# ambiguity refusal and expect_silent tested below ---
muInteractionsSpec <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(forest(interactions = interactions(max.order = 1L)))
)
muInteractionAttr <-
  attr(muInteractionsSpec$control, "bartcore.forests")$interactions[[1L]]
expect_equal(muInteractionAttr$max.order, 1L)

tauInteractionsSpec <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(
    forest(),
    forest(interactions = interactions(max.order = 2L))
  )
)
tauInteractionAttr <-
  attr(tauInteractionsSpec$control, "bartcore.forests")$interactions[[2L]]
expect_equal(tauInteractionAttr$max.order, 2L)

tauBlocksSpec <- dbartsSpec(
  dbartsData(x, y, bases = list(NULL, zBasis)),
  seededControl(),
  forests = list(
    forest(),
    forest(blocks = blocks(groups = list(c("x1", "x2"), c("x3", "x4"))))
  )
)
tauBlockAttr <- attr(tauBlocksSpec$control, "bartcore.forests")$blocks[[2L]]
expect_equal(tauBlockAttr$block.of.column, c(0L, 0L, 1L, 1L))
expect_equal(tauBlockAttr$block.tree.counts, c(25L, 25L))

# --- the creation-refusal contract, inherited half: every option the
# two-forest chain does not read refuses at creation, one assertion each. A
# silently accepted one is the failure mode this creation route exists to
# prevent. ---
expect_error(
  dbarts(x, y, forests = twoForests, tree.prior = dart, control = control),
  "DART tree prior"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    tree.prior = cgm(split.probs = c(0.4, 0.2, 0.2, 0.2)),
    control = control
  ),
  "split.probs"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    monotone = c(1, 0, 0, 0),
    control = control
  ),
  "monotone"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    node.prior = linear(columns = "x1"),
    control = control
  ),
  "linear node prior"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    node.prior = gp(columns = "x1"),
    control = control
  ),
  "Gaussian-process node prior"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    node.prior = normal(chi(1.5)),
    control = control
  ),
  "'k' hyperprior"
)
expect_error(
  dbarts(x, y, forests = twoForests, node.prior = normal(3), control = control),
  "non-default 'k'"
)
expect_error(
  dbarts(x, y, forests = twoForests, variance = TRUE, control = control),
  "'variance'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControl(
      storage = "single"
    )
  ),
  "storage = \"single\""
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0.1,
      change = 0.3,
      birth = 0.5
    ),
    control = control
  ),
  "non-default 'proposal.probs'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    resid.dist = student(5),
    control = control
  ),
  "Student-t residuals"
)
expect_error(
  dbarts(x, y, test = x, forests = twoForests, control = control),
  "test predictors"
)
# per-column cut counts survive only through dbartsSpec(), which keeps a data
# object's own resolved n.cuts
unevenCuts <- dbartsData(x, y, bases = list(NULL, zBasis))
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
  dbartsSpec(dbartsData(x, y, bases = list(NULL, zBasis)), groupedControl),
  "grouped random effects"
)
# the doors: a family whose own parameter block is not shown to interleave with
# the amplitude block refuses by name. (The binary families are no longer among
# them - they build, and test-bcf-family.R is where that surface is pinned.)
expect_error(
  dbarts(
    x,
    rpois(n, 2),
    forests = twoForests,
    family = "nbinom",
    control = control
  ),
  "does not support family \"nbinom\""
)

# --- the creation-refusal contract, new half: every declaration today's
# engine cannot honour refuses at creation, one assertion each, on the same
# refuse-rather-than-drop discipline the inherited set above enforces ---
# --- six of the new half's refusals have since become POSITIVE routes.
# Each is asserted by what it now BUILDS - the forest count, the basis
# widths, and the ragged amplitude vector those widths imply - so a
# relaxation that admitted the declaration while dropping it silently would
# still fail here ---
threeForests <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z)),
    forest(basis = x[, 1L])
  ),
  control = control
)
expect_equal(length(threeForests$data@bases), 3L)
expect_equal(ncol(threeForests$data@bases[[3L]]), 1L)
expect_equal(nrow(threeForests$getForestAmplitudes(3L)), 1L)
expect_equal(dim(threeForests$run(0L, 3L)$glue), c(4L, 3L, 2L))

# a CONTINUOUS basis column is legal: one amplitude, no indicator structure
continuousBasis <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = z)),
  control = control
)
expect_equal(ncol(continuousBasis$data@bases[[2L]]), 1L)
expect_equal(nrow(continuousBasis$getForestAmplitudes(2L)), 1L)

# a factor wider than two levels expands to one amplitude PER LEVEL
threeLevel <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = factor(rep(c("a", "b", "c"), length.out = n)))
  ),
  control = control
)
expect_equal(ncol(threeLevel$data@bases[[2L]]), 3L)
expect_equal(nrow(threeLevel$getForestAmplitudes(2L)), 3L)
expect_true(all(rowSums(threeLevel$data@bases[[2L]]) == 1))
# the FIRST forest takes a basis too, and its amplitude block widens with it
firstForestBasis <- dbarts(
  x,
  y,
  forests = list(forest(basis = ~ factor(z)), forest(basis = ~ factor(z))),
  control = control
)
expect_equal(nrow(firstForestBasis$getForestAmplitudes(1L)), 2L)
expect_equal(dim(firstForestBasis$run(0L, 3L)$glue), c(4L, 3L, 2L))

# 'vars' restricts ANY forest, not only a basis one
firstForestVars <- dbarts(
  x,
  y,
  forests = list(forest(vars = "x1"), forest(basis = ~ factor(z))),
  control = control
)
expect_equal(
  attr(firstForestVars$control, "bartcore.forests")$vars[[1L]],
  1L
)

# 'amplitude.prior.variance' is legal wherever a basis is, the first forest
# now included - and refused, restated per forest, where none is
priorOnFirst <- dbarts(
  x,
  y,
  forests = list(
    forest(basis = ~ factor(z), amplitude.prior.variance = 1),
    forest(basis = ~ factor(z))
  ),
  control = control
)
expect_equal(
  attr(priorOnFirst$control, "bartcore.forests")$params[[1L]][6L],
  1
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(amplitude.prior.variance = 1),
      forest(basis = ~ factor(z))
    ),
    control = control
  ),
  "forest 1 has no 'basis'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(basis = ~ factor(z)),
      forest(basis = ~ factor(z)),
      forest(basis = x[, 1L], amplitude.prior.variance = 2),
      forest(amplitude.prior.variance = 3)
    ),
    control = control
  ),
  "forest 4 has no 'basis'"
)
expect_error(
  dbarts(x, y, forests = list(forest(sd = 1.5)), control = control),
  "single-forest 'forests' has none"
)

# --- K = 1 THROUGH THE K-FOREST PATH. A lone forest CARRYING A BASIS is a
# different declaration from the one just above, and the two creation routes
# used to do different and separately wrong things with it: the data
# route refused in the data object's own vocabulary, and the forests route
# reached no error at all - forestBasisDeclarations' two-forest floor meant the
# basis never reached data@bases, so an ordinary single-forest model was fit
# with the declaration SILENTLY DROPPED. Both now reach one designed refusal at
# the site both pass through, which is what lets it name the RESOLVED count and
# where that count came from. ---
expect_error(
  dbarts(x, y, forests = list(forest(basis = ~ factor(z))), control = control),
  "needs at least two forests"
)
expect_error(
  dbarts(dbartsData(x, y, bases = list(zBasis)), control = control),
  "needs at least two forests"
)
expect_error(
  dbartsSpec(
    dbartsData(x, y),
    control,
    forests = list(forest(basis = ~ factor(z)))
  ),
  "needs at least two forests"
)
# the count's SOURCE is named, because a length-1 declaration over a data object
# already carrying two bases replaces them and resolves to one - telling that
# caller they wrote one basis would be false
expect_error(
  dbartsSpec(
    dbartsData(x, y, bases = list(NULL, zBasis)),
    control,
    forests = list(forest(basis = ~ factor(z)))
  ),
  "'basis' declarations resolve to 1"
)
expect_error(
  dbarts(dbartsData(x, y, bases = list(zBasis)), control = control),
  "data object carries 1"
)
# a length-1 forests carrying BOTH a basis and an amplitude knob now answers the
# K = 1 refusal rather than "a single-forest 'forests' has none": hasBasis is
# TRUE, so the amplitude-knob branch no longer fires. A message change, not an
# acceptance change - the call was refused before and is refused now.
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(basis = ~ factor(z), sd = 1.5)),
    control = control
  ),
  "needs at least two forests"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(update.amplitude = FALSE)),
    control = control
  ),
  "single-forest 'forests' has none"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(interactions = interactions(max.order = 1L))),
    interactions = interactions(max.order = 2L),
    control = control
  ),
  "both at the top level and on the first forest"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(
      blocks = blocks(
        groups = list(
          c("x1", "x2"),
          c(
            "x3",
            "x4"
          )
        )
      )
    )),
    blocks = blocks(groups = list(c("x1", "x2"), c("x3", "x4"))),
    control = control
  ),
  "both at the top level and on the first forest"
)
# a forest past the first with no basis, and no basis reaching the data any
# other way, is not a distinct forest at all - RESTATED per forest under K, so
# the refusal names the one that is missing rather than "the second"
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(n.trees = 25L)),
    control = control
  ),
  "forest 2 needs a 'basis'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z)), forest()),
    control = control
  ),
  "forest 3 needs a 'basis'"
)
# the shape of the argument itself
expect_error(
  dbarts(x, y, forests = forest(), control = control),
  "list of forest\\(\\) specifications"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), list(n.trees = 25L)),
    control = control
  ),
  "list of forest\\(\\) specifications"
)
expect_error(
  dbarts(x, y, forests = list(), control = control),
  "'forests' is empty"
)
# a top-level interactions/blocks with NO forest-level twin still constrains
# the first forest, so the ambiguity refusal is not a blanket one
expect_silent(dbarts(
  x,
  y,
  forests = twoForests,
  interactions = interactions(max.order = 1L),
  control = control
))

# --- forest() validates its own knobs at fit time ---
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), base = 1.5)),
    control = control
  ),
  "forest 'base'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), n.trees = 0L)),
    control = control
  ),
  "forest 'n.trees'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), power = -1)),
    control = control
  ),
  "forest 'power'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), sd = 0)),
    control = control
  ),
  "forest 'sd'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), amplitude.prior.variance = 0)
    ),
    control = control
  ),
  "forest 'amplitude.prior.variance'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), update.amplitude = NA)
    ),
    control = control
  ),
  "forest 'update.amplitude'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = y ~ factor(z))),
    control = control
  ),
  "must be one-sided"
)

# --- a creation the factory would refuse raises an R condition, and no
# usable handle escapes. Driven past the R-layer refusal by hand-attaching the
# variance attribute to an already-resolved BCF spec, so the bridge's own
# backstop is what answers. ---
resolved <- dbartsSpec(dbartsData(x, y), seededControl(), forests = twoForests)
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

# --- the two halves of the specification are cross-checked in BOTH
# directions, so a stripped one is a loud refusal naming the missing piece
# rather than a silent single-forest fit ---
strippedConfig <- resolved
attr(strippedConfig$control, "bartcore.forests") <- NULL
expect_error(
  new(
    "dbartsSampler",
    strippedConfig$control,
    strippedConfig$model,
    strippedConfig$data
  ),
  "no basis forest was configured"
)
strippedTreatment <- resolved
strippedTreatment$data@bases <- NULL
expect_error(
  new(
    "dbartsSampler",
    strippedTreatment$control,
    strippedTreatment$model,
    strippedTreatment$data
  ),
  "carry no forest bases"
)

# --- the pre-built-data refusal: a basis declaration on 'forests' has
# nowhere to ride once 'formula' is already a built dbartsData - dbarts()
# refuses that combination by name rather than silently discarding the
# declaration and fitting a single-forest model; dbartsSpec() is not
# touched, since its first argument is always a dbartsData already; and the
# supported composition - bases on the data object, a knob-only 'forests' -
# keeps working either way. ---
prebuiltData <- dbartsData(x, y)

# (i) the refusal fires by name, on every shape of a basis-declaring
# 'forests' reaching an already-built data object
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = control
  ),
  "cannot reach a pre-built data object"
)
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(basis = ~ factor(z)), forest()),
    control = control
  ),
  "cannot reach a pre-built data object"
)
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(basis = ~ factor(z)), forest(basis = ~ factor(z))),
    control = control
  ),
  "cannot reach a pre-built data object"
)
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(basis = ~ factor(z))),
    control = control
  ),
  "cannot reach a pre-built data object"
)
# the message names both doors out
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = control
  ),
  "dbartsSpec"
)
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = control
  ),
  "dbartsData\\(bases = \\)"
)
# a knob-only 'forests' over a pre-built data object carrying NO bases at all
# is the ordinary single-forest route and is unaffected: no basis is declared
# for the predicate to catch
expect_identical(
  dbarts(prebuiltData, forests = list(forest()), control = control)$run(
    0L,
    3L
  )$train,
  dbarts(prebuiltData, control = control)$run(0L, 3L)$train
)

# (ii) dbartsSpec() takes the SAME declaration, untouched - its first argument
# is always a dbartsData, so the predicate above is unconditionally true there
# and would refuse the one route this surface exists to support
specBuild <- dbartsSpec(
  prebuiltData,
  control,
  forests = list(forest(), forest(basis = ~ factor(z)))
)
specSampler2 <- do.call(
  function(control, model, data, family) {
    new("dbartsSampler", control, model, data)
  },
  specBuild
)
specResult2 <- specSampler2$run(0L, 3L)
expect_equal(dim(specResult2$forestFits), c(n, 2L, 3L, 2L))
expect_equal(dim(specResult2$glue), c(3L, 3L, 2L))
# and the declaration reaches the same model as the ordinary x/y route -
# dbartsSpec() installs it rather than dropping it, bitwise
expect_identical(
  specResult2$train,
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = control
  )$run(0L, 3L)$train
)
# a basis on the first forest too, mirroring (i)'s second refusal, still
# builds through dbartsSpec()
firstForestSpecBuild <- dbartsSpec(
  prebuiltData,
  control,
  forests = list(forest(basis = ~ factor(z)), forest(basis = ~ factor(z)))
)
expect_equal(firstForestSpecBuild$data@bases[[1L]], zBasis)
# the K = 1 route reaches the pre-existing "needs at least two forests"
# refusal here, not the pre-built-data one - a different message, so the two
# refusals are not conflated
expect_error(
  dbartsSpec(
    prebuiltData,
    control,
    forests = list(forest(basis = ~ factor(z)))
  ),
  "needs at least two forests"
)

# (iii) the supported composition - bases already on the data object, a
# knob-only 'forests' - is unaffected: forestBasisDeclarations() returns NULL
# entries for a forest declaring no basis of its own
supportedData <- dbartsData(x, y, bases = list(NULL, zBasis))
supportedResult <- dbarts(
  supportedData,
  forests = list(forest(), forest()),
  control = control
)$run(0L, 3L)
expect_true(!is.null(supportedResult$forestFits))
expect_equal(dim(supportedResult$forestFits), c(n, 2L, 3L, 2L))
expect_equal(dim(supportedResult$glue), c(3L, 3L, 2L))
# and reaches the same model as declaring the basis directly, bitwise
expect_identical(
  supportedResult$train,
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z))),
    control = control
  )$run(0L, 3L)$train
)
# a forests = declaration may still carry KNOBS over that composition - only
# a basis is refused, not the whole argument
supportedKnobSampler <- dbarts(
  supportedData,
  forests = list(forest(), forest(n.trees = 30L)),
  control = control
)
expect_equal(
  attr(supportedKnobSampler$control, "bartcore.forests")$params[[2L]][1L],
  30
)

# (iv) dbartsData()'s own ignored-argument warning now names 'bases' too
warnings.bases <- captureWarnings(
  dbartsData(prebuiltData, bases = list(NULL, zBasis))
)
expect_equal(length(warnings.bases), 1L)
expect_match(conditionMessage(warnings.bases[[1L]]), "ignored")
expect_identical(
  suppressWarnings(dbartsData(prebuiltData, bases = list(NULL, zBasis))),
  prebuiltData
)
expect_silent(dbartsData(prebuiltData))
expect_identical(dbartsData(prebuiltData), prebuiltData)
# the warning is specific to the dbartsData-inherits branch; the ordinary
# x/y route's own 'bases' handling is untouched
expect_silent(dbartsData(x, y, bases = list(NULL, zBasis)))

# (v) the amplitude-prior check now honours 'hasBasis', not only a forest's
# own 'basis' argument: a forest excused by a basis that arrived through
# dbartsData(bases = ) accepts 'amplitude.prior.variance', and one with no
# basis anywhere is still refused
dataBasisOnFirst <- dbartsData(x, y, bases = list(zBasis, zBasis))
amplitudeExcused <- dbarts(
  dataBasisOnFirst,
  forests = list(forest(amplitude.prior.variance = 2), forest()),
  control = control
)
expect_equal(
  attr(amplitudeExcused$control, "bartcore.forests")$params[[1L]][6L],
  2
)
# forest 2's own default applies too, unaffected by forest 1's excuse
expect_equal(
  attr(amplitudeExcused$control, "bartcore.forests")$params[[2L]][6L],
  0.5
)
dataNoBasisOnFirst <- dbartsData(x, y, bases = list(NULL, zBasis))
expect_error(
  dbarts(
    dataNoBasisOnFirst,
    forests = list(forest(amplitude.prior.variance = 2), forest()),
    control = control
  ),
  "forest 1 has no 'basis'"
)
# the excuse applies per forest, not only to the first: forest 2's own
# amplitude.prior.variance is honoured too when ITS basis arrived via the
# data object
amplitudeExcusedBoth <- dbarts(
  dataBasisOnFirst,
  forests = list(forest(), forest(amplitude.prior.variance = 3)),
  control = control
)
expect_equal(
  attr(amplitudeExcusedBoth$control, "bartcore.forests")$params[[2L]][6L],
  3
)
# forest 1's excuse does not extend to a second forest that has no basis
# anywhere either - the exclusion is per forest, not blanket
dataBasisOnFirstOnly <- dbartsData(x, y, bases = list(zBasis, NULL))
expect_error(
  dbarts(
    dataBasisOnFirstOnly,
    forests = list(forest(amplitude.prior.variance = 2), forest()),
    control = control
  ),
  "forest 2 needs a 'basis'"
)
# and an amplitude.prior.variance declared on THAT unexcused forest answers
# the amplitude refusal instead, since the amplitude-prior check runs first
expect_error(
  dbarts(
    dataBasisOnFirstOnly,
    forests = list(forest(), forest(amplitude.prior.variance = 2)),
    control = control
  ),
  "forest 2 has no 'basis'"
)
