# Public multi-forest creation: dbarts() and dbartsSpec()
# take a forests = list(forest(...)) declaration and build an ordinary
# dbartsSampler holding the model it names. Two forests, the second carrying a
# two-level factor basis, is the Bayesian causal forest. Every option
# the two-forest chain does not read, and every declaration today's engine
# cannot honor, must refuse at creation rather than be dropped in silence.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

source(
  system.file("common", "forestParams.R", package = "dbarts"),
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

uncountedControlBcfCreation <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    updateState = FALSE,
    seed = 17L,
    ...
  )
}
seededControlBcfCreation <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    updateState = FALSE,
    seed = 17L,
    ...
  )
}

# the declaration every refusal below is attached to: a plain first forest and
# a second one whose two-level factor basis carries the (b0, b1) amplitudes
twoForests <- list(
  dbarts::dbartsForests$forest(),
  dbarts::dbartsForests$forest(basis = ~ factor(z))
)

# --- creation: with control@seed set, every chain's generator is independent
# of R's stream. The removed treatment = route needs no reconstruction. ---
control <- seededControlBcfCreation()
publicSampler <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z), n.trees = 25L, sd = 1.5)
  ),
  control = control,
  n.trees = 50L
)
publicResult <- publicSampler$run(0L, 10L)

# the sampler this produced is an ordinary dbartsSampler carrying two forests,
# and the level indicators the factor basis expanded to are the ones the slot
# takes, forest by forest - the first forest declaring none
expect_true(inherits(publicSampler, "dbartsSampler"))
expect_null(publicSampler$data@bases[[1L]])
expect_equal(
  publicSampler$data@bases[[2L]],
  cbind("factor(z)0" = 1 - as.double(z), "factor(z)1" = as.double(z))
)
expect_true(all(is.finite(publicResult$train)))
expect_true(all(publicResult$sigma > 0))

# --- the expansion's level ORDER is load-bearing. Swapping the two indicator
# columns swaps which rows b1 scales, so a sampler built from the reversed
# factor must NOT reproduce the un-swapped one. ---
swappedSampler <- dbarts(
  x,
  y,
  forests = list(
    forest(),
    forest(basis = ~ factor(z, levels = c(1, 0)), n.trees = 25L, sd = 1.5)
  ),
  control = seededControlBcfCreation(),
  n.trees = 50L
)
expect_equal(
  swappedSampler$data@bases[[2L]],
  cbind(
    "factor(z, levels = c(1, 0))1" = as.double(z),
    "factor(z, levels = c(1, 0))0" = 1 - as.double(z)
  )
)
expect_false(identical(
  swappedSampler$run(0L, 10L)$train,
  publicResult$train
))

# --- the reported draws really are the two-forest blend: the amplitudes and
# the per-forest fits reconstruct the recorded train draw through the stored
# response transform (the identity a low-level pin already established, now
# on a public sampler) ---
glue <- publicSampler$getForestAmplitudes()
muFits <- publicSampler$getForestFits(1L)[, 1L]
tauFits <- publicSampler$getForestFits(2L)[, 1L]
publicSampler$storeState()
fitScale <- publicSampler$state[[1L]]$fit.scale
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
  control = seededControlBcfCreation(),
  n.trees = 50L
)
restricted$run(0L, 10L)
tauCounts <- restricted$getForestVariableCounts(2L)
expect_true(all(tauCounts[1:2, ] == 0L))
expect_true(sum(tauCounts[3:4, ]) > 0L)

# --- the whole knob map is honored, not accepted and dropped. Each knob of
# the second forest reaches the four doubles the control attribute carries
# (its `sd` goes beside them, as written), and
# the FIRST forest's structural knobs restate the fit's own control and tree
# prior rather than adding a second set. A knob that went nowhere would leave
# a default here. ---
knobs <- dbartsSpec(
  dbartsData(x, y),
  uncountedControlBcfCreation(),
  forests = list(
    forest(
      n.trees = 30L,
      base = 0.7,
      power = 1.5,
      sd = 2.5,
      amplitude = fixed()
    ),
    forest(
      basis = ~ factor(z),
      n.trees = 15L,
      base = 0.4,
      power = 2,
      sd = 1.25,
      amplitude = fixed()
    )
  )
)
# one length-4 vector per forest: tree count, base, power and the hold, 0 where
# the coefficient is held. Which magnitude channel a forest's `sd` reaches is
# the engine's to decide from whether it carries a basis, so the stated sds go
# to the engine as written, in the response's units, one a forest; no forest
# states a variance for its coefficients
expect_equal(
  forestParamsOf(knobs$model),
  list(c(30, 0.7, 1.5, 0), c(15, 0.4, 2, 0))
)
expect_identical(vapply(knobs$model@forests, slot, 0, "sd"), c(2.5, 1.25))
expect_identical(knobs$control@n.trees, 30L)
expect_equal(knobs$model@tree.prior@base, 0.7)
expect_equal(knobs$model@tree.prior@power, 1.5)

# --- a forest() n.trees on the first forest is a more specific declaration
# than an explicitly supplied control@n.trees, and governs the fit over it
# rather than being silently overridden by it: the run really is the
# seven-tree fit, not the hundred-tree one it disagrees with. Its base and
# power are the fitting function's own arguments, which beside it are one
# setting written twice and refused ---
sevenTreeControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 7L, base = 0.3, power = 1.1)),
    control = uncountedControlBcfCreation()
  )$run(0L, 5L)$train,
  dbarts(
    x,
    y,
    base = 0.3,
    power = 1.1,
    control = sevenTreeControl,
    n.trees = 7L
  )$run(0L, 5L)$train
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 7L, base = 0.3, power = 1.1)),
    base = 0.6,
    control = uncountedControlBcfCreation()
  ),
  paste0(
    "'base' is given to the fitting function and to the forest with no ",
    "basis (forest 1), which are the same setting; give one"
  ),
  fixed = TRUE
)

# The pinned-amplitude public fit: amplitude = fixed() on both forests is
# the pinned-amplitude model (y = mu + z tau exactly) the on-ramp vignette's
# continuity falsifier composes against; its glue is exactly (1, 0, 1)
# (test-bcf.R's bcFixed pins the value, this just confirms the public forest()
# knob reaches the same fixed-glue path)
pinnedPublic <- dbarts(
  x,
  y,
  forests = list(
    forest(amplitude = fixed()),
    forest(basis = ~ factor(z), n.trees = 25L, amplitude = fixed())
  ),
  control = seededControlBcfCreation(),
  n.trees = 50L
)
pinnedPublic$run(0L, 10L)
expect_equal(pinnedPublic$getForestAmplitudes()[, 1L], c(1, 0, 1))

# --- forests = NULL is byte-neutral, and a single-forest declaration is
# the same fit with its structural knobs restated ---
plainResult <- dbarts(
  x,
  y,
  control = seededControlBcfCreation(),
  n.trees = 50L
)$run(0L, 5L)
expect_identical(
  dbarts(
    x,
    y,
    forests = NULL,
    control = seededControlBcfCreation(),
    n.trees = 50L
  )$run(
    0L,
    5L
  )$train,
  plainResult$train
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 50L)),
    control = uncountedControlBcfCreation()
  )$run(0L, 5L)$train,
  plainResult$train
)
# and a restated knob really is the fit's own: a different tree count draws
# what the control carrying that count draws
twentyTrees <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L
)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(n.trees = 20L)),
    control = uncountedControlBcfCreation()
  )$run(0L, 5L)$train,
  dbarts(x, y, control = twentyTrees, n.trees = 20L)$run(0L, 5L)$train
)

# --- dbartsSpec() reaches the same model, and carries the basis column on the
# data object rather than on the control ---
specSampler <- do.call(
  # the spec's own elements by name, the trailing ones (family, the active-row
  # mask a latent family's 0/1 weights resolve to) absorbed rather than listed
  function(control, model, data, ...) {
    new("dbartsSampler", control, model, data)
  },
  dbartsSpec(
    dbartsData(x, y),
    seededControlBcfCreation(),
    forests = list(
      forest(),
      forest(basis = ~ factor(z), n.trees = 25L, sd = 1.5)
    ),
    n.trees = 50L
  )
)
expect_identical(specSampler$run(0L, 10L)$train, publicResult$train)

# a basis handed to a forest beside a data object is restricted by 'subset'
# exactly as the weights beside it are
zFactor <- factor(z)
# the block the factor is stored as, its indicator columns
zBasis <- cbind(1 - as.double(z), as.double(z))
# a sampler over a data object whose second forest multiplies a value
valueBesideData <- function(value, ...) {
  dbarts(
    dbartsData(x, y, ...),
    forests = list(forest(), do.call(forest, list(basis = value))),
    control = seededControlBcfCreation(),
    n.trees = 50L
  )
}
subsetData <- valueBesideData(zFactor, subset = 1:100)$data
expect_equal(unname(subsetData@bases[[2L]]), zBasis[1:100, ])
expect_error(
  valueBesideData(cbind(NA_real_, 1)),
  "forest 2's 'basis' has 1 row where the data were handed",
  fixed = TRUE
)
expect_error(
  valueBesideData(zFactor[1:10]),
  "forest 2's 'basis' has 10 rows where the data were handed",
  fixed = TRUE
)
# a non-finite value is refused rather than reaching the multiplier
expect_error(
  valueBesideData(rep(Inf, n)),
  "a 'basis' must be finite",
  fixed = TRUE
)

# ... while a caller who declared a forest's BASIS is refused in that word, on
# every surface that takes one, rather than in the data object's
shortBasisRefusal <- paste0(
  "'basis' (factor(z[1:10])) must have the same length as the data: it has ",
  "10 rows and the data ",
  length(y)
)
expect_error(
  dbartsSpec(
    dbartsData(x, y),
    seededControlBcfCreation(),
    forests = list(forest(), forest(basis = ~ factor(z[1:10]))),
    n.trees = 50L
  ),
  shortBasisRefusal,
  fixed = TRUE
)
expect_error(
  dbarts(x, y, forests = list(forest(), forest(basis = ~ factor(z[1:10])))),
  shortBasisRefusal,
  fixed = TRUE
)
# a value handed over is refused in the same word
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), do.call(forest, list(basis = factor(z[1:10]))))
  ),
  "forest 2's 'basis' has 10 rows where the data were handed",
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    control = seededControlBcfCreation(),
    forests = list(forest(), forest(basis = ~ factor(z))),
    n.trees = 50L
  )$setForestBasis(2L, ~ factor(z[1:10])),
  shortBasisRefusal,
  fixed = TRUE
)
offSlot <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(forest(), forest(basis = zFactor)),
  n.trees = 50L
)
# the transport is RAGGED - one length-4 vector per forest, not one length-4
# vector for the pair - and a forest that states no sd or variance sends NA for
# each, the engine holding every default
expect_equal(
  forestParamsOf(offSlot$model),
  list(c(50, 0.95, 2, 1), c(50, 0.25, 3, 1))
)
expect_identical(
  vapply(offSlot$model@forests, slot, 0, "sd"),
  c(NA_real_, NA_real_)
)
# the same transport under a LATENT family, which no assertion above has
# run: each forest's `sd` is carried in its own position as written, and the
# family reaches the model without disturbing either. The two declared values
# DIFFER, and neither is 1, so a confusion between a forest's own entry and its
# neighbor's is visible.
latentSlots <- dbartsSpec(
  dbartsData(x, as.double(y > median(y))),
  seededControlBcfCreation(),
  forests = list(forest(sd = 3.5), forest(basis = zFactor, sd = 1.25)),
  family = "logistic",
  n.trees = 50L
)
expect_equal(latentSlots$model@family, "logistic")
expect_equal(
  forestParamsOf(latentSlots$model),
  list(c(50, 0.95, 2, 1), c(50, 0.25, 3, 1))
)
expect_identical(vapply(latentSlots$model@forests, slot, 0, "sd"), c(3.5, 1.25))
# and a basis declared alongside subset reaches dbartsData's own alignment
subsetForests <- dbarts(
  x,
  y,
  subset = 1:100,
  forests = twoForests,
  control = seededControlBcfCreation(),
  n.trees = 50L
)
expect_equal(
  subsetForests$data@bases[[2L]],
  structure(
    zBasis[1:100, ],
    dimnames = list(NULL, c("factor(z)0", "factor(z)1"))
  )
)

# --- the knobs of a forest with no basis beside a forest with one really
# reach the resolved spec ---
sdViaHasBasis <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(forest(sd = 3.5), forest(basis = zFactor)),
  n.trees = 50L
)
expect_identical(vapply(sdViaHasBasis$model@forests, slot, 0, "sd")[1L], 3.5)

noUpdateViaHasBasis <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(forest(amplitude = fixed()), forest(basis = zFactor)),
  n.trees = 50L
)
expect_equal(
  forestParamsOf(noUpdateViaHasBasis$model)[[1L]][4L],
  0
)

# a second forest with no basis takes its own knobs
noBasisSecondForest <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(forest(), forest(n.trees = 25L)),
  n.trees = 50L
)
expect_equal(
  forestParamsOf(noBasisSecondForest$model)[[2L]][1L],
  25
)

# --- interactions()/blocks() declared on a forest() reach the resolved
# per-forest constraint the control attribute carries, not only the top-level
# ambiguity refusal and expect_silent tested below ---
muInteractionsSpec <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(
    forest(interactions = interactions(max.order = 1L)),
    forest(basis = zFactor)
  ),
  n.trees = 50L
)
muInteractionAttr <-
  muInteractionsSpec$model@forests[[1L]]@interactions
expect_equal(muInteractionAttr$max.order, 1L)

tauInteractionsSpec <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(
    forest(),
    forest(basis = zFactor, interactions = interactions(max.order = 2L))
  ),
  n.trees = 50L
)
tauInteractionAttr <-
  tauInteractionsSpec$model@forests[[2L]]@interactions
expect_equal(tauInteractionAttr$max.order, 2L)

tauBlocksSpec <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = list(
    forest(),
    forest(
      basis = zFactor,
      blocks = blocks(groups = list(c("x1", "x2"), c("x3", "x4")))
    )
  ),
  n.trees = 50L
)
tauBlockAttr <- tauBlocksSpec$model@forests[[2L]]@blocks
expect_equal(tauBlockAttr$block.of.column, c(0L, 0L, 1L, 1L))
expect_equal(tauBlockAttr$block.tree.counts, c(25L, 25L))

# --- the creation-refusal contract, inherited half: every option the
# two-forest chain does not read refuses at creation, one assertion each. A
# silently accepted one is the failure mode this creation route exists to
# prevent. ---
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    split.probs = dart,
    control = control,
    n.trees = 50L
  ),
  "does not support 'split.probs' (a vector that is not uniform, or dart())",
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    split.probs = c(0.4, 0.2, 0.2, 0.2),
    control = control,
    n.trees = 50L
  ),
  "split.probs"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    monotone = c(1, 0, 0, 0),
    control = control,
    n.trees = 50L
  ),
  "monotone"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    leaf.prior = linear(columns = "x1"),
    control = control,
    n.trees = 50L
  ),
  "linear leaf prior"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    leaf.prior = gp(columns = "x1"),
    control = control,
    n.trees = 50L
  ),
  "Gaussian-process leaf prior"
)
# with a forest under it the fitting function's leaf prior is taken; with none
# it has no forest
noneTakes <- list(
  dbartsForests$forest(sd = 1),
  dbartsForests$forest(basis = ~ factor(z))
)
expect_error(
  dbarts(
    x,
    y,
    forests = noneTakes,
    leaf.prior = normal(chi(1.5)),
    control = control,
    n.trees = 50L
  ),
  "takes a leaf prior's 'k'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = noneTakes,
    leaf.prior = normal(3),
    control = control,
    n.trees = 50L
  ),
  "takes a leaf prior's 'k'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    variance = TRUE,
    control = control,
    n.trees = 50L
  ),
  "'variance'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControlBcfCreation(
      storage = "single"
    ),
    n.trees = 50L
  ),
  "storage = \"single\""
)
# the tree-move mixture is NOT on the refusal list: it is a property of the
# fit, carried onto every forest of the declaration, so each of these creates.
# That the mixture reaches BOTH forests rather than only the prognostic one is
# pinned in test-proposal-probs.R, which reads the two forests' split counts.
expect_silent(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControlBcfCreation(
      proposal.probs = c(
        birth_death = 0.6,
        swap = 0.1,
        change = 0.3,
        birth = 0.5
      )
    ),
    n.trees = 50L
  )
)
expect_silent(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControlBcfCreation(
      proposal.probs = c(
        birth_death = 0.6,
        swap = 0,
        change = 0.4,
        birth = 0.5
      )
    ),
    n.trees = 50L
  )
)
expect_silent(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControlBcfCreation(
      proposal.probs = c(
        birth_death = 0.6,
        swap = 0,
        change = 0.24,
        perturb = 0.16,
        birth = 0.5
      )
    ),
    n.trees = 50L
  )
)
expect_silent(
  dbarts(
    x,
    y,
    forests = twoForests,
    control = seededControlBcfCreation(
      proposal.probs = c(
        birth_death = 0.6,
        swap = 0,
        change = 0.24,
        rule_gibbs = 0.16,
        birth = 0.5
      )
    ),
    n.trees = 50L
  )
)
expect_error(
  dbarts(
    x,
    y,
    forests = twoForests,
    family = student(5),
    control = control,
    n.trees = 50L
  ),
  "Student-t residuals"
)
expect_error(
  dbarts(
    x,
    y,
    test = x,
    forests = twoForests,
    control = control,
    n.trees = 50L
  ),
  "test predictors"
)
# per-column cut counts survive only through dbartsSpec(), which keeps a data
# object's own resolved n.cuts
unevenCuts <- dbartsData(x, y, n.cuts = c(50L, 100L, 100L, 100L))
expect_error(
  dbartsSpec(
    unevenCuts,
    seededControlBcfCreation(),
    forests = list(forest(), forest(basis = zFactor)),
    n.trees = 50L
  ),
  "per-column 'n.cuts'"
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
    control = control,
    n.trees = 50L
  ),
  "does not support family \"nbinom\""
)

# --- the creation-refusal contract, new half: every declaration today's
# engine cannot honor refuses at creation, one assertion each, on the same
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
  control = control,
  n.trees = 50L
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
  control = control,
  n.trees = 50L
)
expect_equal(ncol(continuousBasis$data@bases[[2L]]), 1L)
expect_equal(nrow(continuousBasis$getForestAmplitudes(2L)), 1L)

# a factor wider than two levels is refused, a forest having two levels or
# one numeric column
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = factor(rep(c("a", "b", "c"), length.out = n)))
    ),
    control = control,
    n.trees = 50L
  ),
  "has a basis of 3 levels ('a', 'b', 'c')",
  fixed = TRUE
)
# the FIRST forest takes a basis too, and its amplitude block widens with it.
# No forest of the model is without a basis, so the count and the tree prior
# are stated on a forest and the control names no count
firstForestBasis <- dbarts(
  x,
  y,
  forests = list(
    forest(basis = ~ factor(z), n.trees = 50L, base = 0.95, power = 2),
    forest(basis = ~ factor(z))
  ),
  control = uncountedControlBcfCreation()
)
expect_equal(nrow(firstForestBasis$getForestAmplitudes(1L)), 2L)
expect_equal(dim(firstForestBasis$run(0L, 3L)$glue), c(4L, 3L, 2L))

# 'vars' restricts ANY forest, not only a basis one
firstForestVars <- dbarts(
  x,
  y,
  forests = list(forest(vars = "x1"), forest(basis = ~ factor(z))),
  control = control,
  n.trees = 50L
)
expect_equal(
  firstForestVars$model@forests[[1L]]@vars,
  1L
)
# a fractional 'vars' entry is refused, naming the argument, rather than
# silently truncated (coerceOrError's integer branch); a whole-number double
# entry is accepted the same as an integer or name
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(vars = 2.9), forest(basis = ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "'vars' must be a whole number; got '2.9'",
  fixed = TRUE
)
firstForestVarsNumeric <- dbarts(
  x,
  y,
  forests = list(forest(vars = 2), forest(basis = ~ factor(z))),
  control = control,
  n.trees = 50L
)
expect_equal(
  firstForestVarsNumeric$model@forests[[1L]]@vars,
  2L
)

# 'amplitude.prior.variance' is not an argument of forest(): R's own error,
# pinned in test-forest-held.R
# the one forest of a single-forest model takes `sd` as the fitting function's
# leaf.prior = normal(sd = ) does, in the same units: the same fit, draw for
# draw
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(sd = 1.5)),
    control = seededControlBcfCreation(),
    n.trees = 50L
  )$run(0L, 10L),
  dbarts(
    x,
    y,
    leaf.prior = normal(sd = 1.5),
    control = seededControlBcfCreation(),
    n.trees = 50L
  )$run(0L, 10L)
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
  dbarts(
    x,
    y,
    forests = list(forest(basis = ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "needs at least two forests"
)
expect_error(
  dbarts(
    dbartsData(x, y),
    forests = list(do.call(forest, list(basis = zFactor))),
    control = control,
    n.trees = 50L
  ),
  "needs at least two forests"
)
expect_error(
  dbartsSpec(
    dbartsData(x, y),
    control,
    forests = list(forest(basis = ~ factor(z))),
    n.trees = 50L
  ),
  "needs at least two forests"
)
# a length-1 forests carrying BOTH a basis and an amplitude knob now answers the
# K = 1 refusal rather than the one-forest refusal of the amplitude knobs:
# hasBasis is TRUE, so that branch no longer fires. A message change, not an
# acceptance change - the call was refused before and is refused now.
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(basis = ~ factor(z), sd = 1.5)),
    control = control,
    n.trees = 50L
  ),
  "needs at least two forests"
)
# fixed() on a lone forest restates its default, with no message (dec-B513)
set.seed(44L)
loneHeld <- dbarts(
  x,
  y,
  forests = list(forest(amplitude = fixed())),
  control = control,
  n.trees = 50L
)
set.seed(44L)
lonePlain <- dbarts(x, y, control = control, n.trees = 50L)
expect_identical(loneHeld$run(0L, 3L), lonePlain$run(0L, 3L))
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(interactions = interactions(max.order = 1L))),
    interactions = interactions(max.order = 2L),
    control = control,
    n.trees = 50L
  ),
  "is given to the fitting function and to the forest with no basis"
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
    control = control,
    n.trees = 50L
  ),
  "is given to the fitting function and to the forest with no basis"
)
# two forests with no basis share one prior total, wherever they stand
sharedPositions <- function(sampler) {
  which(vapply(
    unname(sampler$getLeafPrior()),
    function(prior) identical(prior$prior.sd.of, "forest, shared"),
    NA
  ))
}
expect_identical(
  sharedPositions(dbarts(
    x,
    y,
    forests = list(forest(), forest(n.trees = 25L)),
    control = control,
    n.trees = 50L
  )),
  1:2
)
expect_identical(
  sharedPositions(dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z)), forest()),
    control = control,
    n.trees = 50L
  )),
  c(1L, 3L)
)
# the shape of the argument itself
expect_error(
  dbarts(x, y, forests = 3L, control = control, n.trees = 50L),
  "list of forest\\(\\) specifications"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), list(n.trees = 25L)),
    control = control,
    n.trees = 50L
  ),
  "list of forest\\(\\) specifications"
)
expect_error(
  dbarts(x, y, forests = list(), control = control, n.trees = 50L),
  "'forests' is empty"
)
# a top-level interactions/blocks with NO forest-level twin still constrains
# the first forest, so the ambiguity refusal is not a blanket one
expect_silent(dbarts(
  x,
  y,
  forests = twoForests,
  interactions = interactions(max.order = 1L),
  control = control,
  n.trees = 50L
))

# --- forest() validates its own knobs at fit time ---
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), base = 1.5)),
    control = control,
    n.trees = 50L
  ),
  "'base'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), n.trees = 0L)),
    control = control,
    n.trees = 50L
  ),
  "'n.trees'"
)
# a fractional n.trees is refused, naming the argument, rather than silently
# truncated (coerceOrError's integer branch)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), n.trees = 5.5)),
    control = control,
    n.trees = 50L
  ),
  "'n.trees' must be a whole number; got '5.5'",
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), power = -1)),
    control = control,
    n.trees = 50L
  ),
  "'power'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~ factor(z), sd = 0)),
    control = control,
    n.trees = 50L
  ),
  "'sd'"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), amplitude = FALSE)
    ),
    control = control,
    n.trees = 50L
  ),
  "'amplitude' must be fixed()"
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = y ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "must be one-sided"
)

# --- a creation the factory would refuse raises an R condition, and no
# usable handle escapes. Driven past the R-layer refusal by hand-attaching the
# variance attribute to an already-resolved BCF spec, so the bridge's own
# backstop is what answers. ---
resolved <- dbartsSpec(
  dbartsData(x, y),
  seededControlBcfCreation(),
  forests = twoForests,
  n.trees = 50L
)
withVariance <- resolved
withVariance$model@variance <- methods::new(
  "dbartsVarianceRecord",
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
scaledModel$model@leaf.scale <- 3.0
expect_error(
  new(
    "dbartsSampler",
    scaledModel$control,
    scaledModel$model,
    scaledModel$data
  ),
  "non-default leaf scale"
)

# --- the two halves of the specification are cross-checked in BOTH
# directions, so a stripped one is a loud refusal naming the missing piece
# rather than a silent single-forest fit ---
strippedConfig <- resolved
strippedConfig$model@forests <- strippedConfig$model@forests[1L]
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
strippedTreatment$data@basis.levels <- NULL
expect_error(
  new(
    "dbartsSampler",
    strippedTreatment$control,
    strippedTreatment$model,
    strippedTreatment$data
  ),
  "carry no forest bases"
)

# --- a pre-built data object: a basis declared on 'forests' is read beside
# it, as at the data it was built from, at dbarts() and dbartsSpec() alike,
# and a list of one forest with a basis is refused as anywhere ---
prebuiltData <- dbartsData(x, y)

# (i) the declaration reaches the same model as the ordinary x/y route, in
# every shape
for (declared in list(
  quote(list(forest(), forest(basis = ~ factor(z)))),
  quote(list(
    forest(basis = ~ factor(z), n.trees = 50L, base = 0.95, power = 2),
    forest()
  )),
  quote(list(
    forest(basis = ~ factor(z), n.trees = 50L, base = 0.95, power = 2),
    forest(basis = ~ factor(z))
  ))
)) {
  expect_identical(
    eval(bquote(dbarts(
      prebuiltData,
      forests = .(declared),
      control = uncountedControlBcfCreation()
    )))$run(0L, 3L)$train,
    eval(bquote(dbarts(
      x,
      y,
      forests = .(declared),
      control = uncountedControlBcfCreation()
    )))$run(0L, 3L)$train,
    info = deparse(declared)
  )
}
expect_error(
  dbarts(
    prebuiltData,
    forests = list(forest(basis = ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "needs at least two forests"
)
# a knob-only 'forests' over a pre-built data object carrying NO bases at all
# is the ordinary single-forest route and is unaffected: no basis is declared
# for the predicate to catch
expect_identical(
  dbarts(
    prebuiltData,
    forests = list(forest()),
    control = control,
    n.trees = 50L
  )$run(
    0L,
    3L
  )$train,
  dbarts(prebuiltData, control = control, n.trees = 50L)$run(0L, 3L)$train
)

# (ii) dbartsSpec() takes the SAME declaration, untouched - its first argument
# is always a dbartsData, so the predicate above is unconditionally true there
# and would refuse the one route this surface exists to support
specBuild <- dbartsSpec(
  prebuiltData,
  control,
  forests = list(forest(), forest(basis = ~ factor(z))),
  n.trees = 50L
)
specSampler2 <- do.call(
  # the spec's own elements by name, the trailing ones (family, the active-row
  # mask a latent family's 0/1 weights resolve to) absorbed rather than listed
  function(control, model, data, ...) {
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
    control = control,
    n.trees = 50L
  )$run(0L, 3L)$train
)
# a basis on the first forest too, mirroring (i)'s second refusal, still
# builds through dbartsSpec()
firstForestSpecBuild <- dbartsSpec(
  prebuiltData,
  uncountedControlBcfCreation(),
  forests = list(
    forest(basis = ~ factor(z), n.trees = 50L, base = 0.95, power = 2),
    forest(basis = ~ factor(z))
  )
)
expect_equal(
  firstForestSpecBuild$data@bases[[1L]],
  structure(zBasis, dimnames = list(NULL, c("factor(z)0", "factor(z)1")))
)
# the K = 1 route reaches the pre-existing "needs at least two forests"
# refusal here, not the pre-built-data one - a different message, so the two
# refusals are not conflated
expect_error(
  dbartsSpec(
    prebuiltData,
    control,
    forests = list(forest(basis = ~ factor(z))),
    n.trees = 50L
  ),
  "needs at least two forests"
)

# (iii) a basis handed over as a value beside the data object
supportedData <- dbartsData(x, y)
supportedResult <- dbarts(
  supportedData,
  forests = list(forest(), forest(basis = zFactor)),
  control = control,
  n.trees = 50L
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
    control = control,
    n.trees = 50L
  )$run(0L, 3L)$train
)
# and with knobs beside it
supportedKnobSampler <- dbarts(
  supportedData,
  forests = list(forest(), forest(basis = zFactor, n.trees = 30L)),
  control = control,
  n.trees = 50L
)
expect_equal(
  forestParamsOf(supportedKnobSampler$model)[[2L]][1L],
  30
)

# (iv) a data object handed to dbartsData() is returned as it is
expect_silent(dbartsData(prebuiltData))
expect_identical(dbartsData(prebuiltData), prebuiltData)

# (v) a forest whose basis is handed over beside a data object is one with a
# basis: its channel and its tree defaults are a multiplied forest's
bothExcused <- dbarts(
  dbartsData(x, y),
  forests = list(
    forest(basis = zFactor, n.trees = 50L, base = 0.95, power = 2),
    forest(basis = zFactor)
  ),
  control = uncountedControlBcfCreation()
)
expect_identical(
  unname(vapply(bothExcused$getLeafPrior(), `[[`, "", "prior.sd.of")),
  rep("level difference (at a coefficient of 1)", 2L)
)
# a second forest that has no basis anywhere is the model's forest with no
# basis, which may stand at any place: it takes the control's count and the
# fitting function's tree prior, and the first forest the multiplied forest's
plainSecond <- dbarts(
  dbartsData(x, y),
  forests = list(forest(basis = zFactor), forest()),
  control = uncountedControlBcfCreation(),
  n.trees = 20L
)
expect_equal(
  forestParamsOf(plainSecond$model),
  list(c(50, 0.25, 3, 1), c(20, 0.95, 2, 1))
)
# each forest's channel follows its kind, not its place: the first, a
# factor, carries its sd in its leaf scale and the second the median
plainSecondPriors <- plainSecond$getLeafPrior()
expect_null(plainSecondPriors[[1L]]$amplitude.prior.scale)
expect_identical(
  plainSecondPriors[[1L]]$prior.sd.of,
  "level difference (at a coefficient of 1)"
)
expect_equal(
  plainSecondPriors[[1L]]$k.scale,
  plainSecondPriors[[1L]]$sd,
  tolerance = 1e-12
)
expect_equal(
  plainSecondPriors[[2L]]$sd,
  plainSecondPriors[[2L]]$k.scale / 2,
  tolerance = 1e-12
)
expect_identical(plainSecondPriors[[2L]]$prior.sd.of, "forest")

# --- dec-A109: the lone basis forest stays refused, and the refusal points at
# the p + 1 forest form and linear() leaves ---
expect_error(
  dbarts(
    dbartsData(x, y),
    forests = list(do.call(forest, list(basis = zFactor))),
    control = control,
    n.trees = 50L
  ),
  "forest(basis = z1)",
  fixed = TRUE
)
expect_error(
  dbarts(
    dbartsData(x, y),
    forests = list(do.call(forest, list(basis = zFactor))),
    control = control,
    n.trees = 50L
  ),
  "linear() leaves",
  fixed = TRUE
)

# --- two basis-free forests beside a data object share one total as the
# forests = route's do ---
expect_identical(
  sharedPositions(dbarts(
    dbartsData(x, y),
    forests = list(forest(), forest()),
    control = control,
    n.trees = 50L
  )),
  1:2
)
spec <- dbartsSpec(
  dbartsData(x, y),
  control,
  forests = list(forest(), forest(), forest(basis = zFactor)),
  n.trees = 50L
)
expect_identical(
  vapply(spec$model@forests, function(record) as.numeric(record@n.trees), 0),
  c(25, 25, 50)
)

# --- xbart and rbart_vi fit one forest: beside a sampler's data, as at
# bart(), no forests list gives one forest and the bases they hold are not part
# of the model ---
basesData <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = zFactor)),
  control = control,
  n.trees = 50L
)$data
crossValidated <- function(data) {
  xbart(
    data,
    n.samples = 5L,
    n.test = 2L,
    n.reps = 2L,
    n.threads = 1L,
    n.trees = 5L,
    seed = 3L
  )
}
expect_identical(crossValidated(basesData), crossValidated(dbartsData(x, y)))
grouped <- function(data) {
  suppressWarnings(rbart_vi(
    data,
    group.by = rep(1:2, length.out = n),
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    seed = 3L,
    verbose = FALSE
  ))$yhat.train
}
expect_identical(grouped(basesData), grouped(dbartsData(x, y)))

# a basis on the FIRST forest and none on the second: the control's slot
# holds the first forest's count and the control's own goes to the second
plainSecondSpec <- dbartsSpec(
  dbartsData(x, y),
  uncountedControlBcfCreation(),
  forests = list(forest(basis = zFactor), forest()),
  n.trees = 20L
)
expect_identical(plainSecondSpec$control@n.trees, 20L)
expect_identical(plainSecondSpec$model@forests[[1L]]@n.trees, 50L)
expect_equal(
  forestParamsOf(plainSecondSpec$model)[[2L]][1:3],
  c(20, 0.95, 2)
)
# the clause names the declaration when the caller wrote forests =
expect_error(
  dbarts(
    y ~ x,
    data.frame(y = y, x = x[, 1L], z = z),
    forests = list(forest(basis = ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "'basis' declarations resolve to 1"
)
expect_error(
  dbarts(
    y ~ x,
    data.frame(y = y, x = x[, 1L], z = z),
    forests = list(forest(basis = ~ factor(z))),
    control = control,
    n.trees = 50L
  ),
  "Write the forest with no multiplier too, as y ~ forest(x1 + x2) + ",
  fixed = TRUE
)
