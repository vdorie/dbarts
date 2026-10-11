# A proposal.probs or n.perturb.cuts changed through setControl reaches both
# forests of a heteroscedastic sampler. The variance forest takes the control's
# mixture and width when the sampler is created, so it takes the new ones with
# the mean forest.

set.seed(31L)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] - x[, 3L] + ifelse(x[, 1L] > 0.5, 0.2, 1) * rnorm(n)
mixtureControl <- function(...) {
  dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE,
    ...
  )
}
heteroscedastic <- function(control, seed) {
  set.seed(seed)
  dbarts::dbarts(
    x,
    y,
    variance = varianceForest(n.trees = 5L),
    control = control
  )
}
# what a run reports and the surface it leaves, the two forests' draws
drawsOf <- function(sampler, seed, n.burn, n.samples) {
  set.seed(seed)
  list(
    train = sampler$run(n.burn, n.samples)$train,
    variance = sampler$getVariance()
  )
}
# variables, cut points and node counts, but not the leaf values
meanStructure <- function(sampler) {
  trees <- sampler$getTrees()
  paste(
    trees$tree,
    trees$n,
    trees$var,
    ifelse(trees$var != -1L, trees$value, NA_real_),
    collapse = "|"
  )
}
# the variance trees' shapes, split variables and cuts, off the stored state:
# a node's value is its cut where it splits and its factor where it does not
varianceStructure <- function(sampler) {
  sampler$storeState()
  chainState <- sampler$state[[1L]]
  splits <- chainState$variance.vars >= 0L
  list(
    vars = chainState$variance.vars,
    sizes = chainState$variance.sizes,
    cuts = matrix(chainState$variance.values, 8L)[, splits, drop = FALSE]
  )
}
sweepsMoved <- function(sampler, structureOf, sweeps = 20L) {
  at <- structureOf(sampler)
  moved <- 0L
  for (sweep in seq_len(sweeps)) {
    invisible(sampler$run(0L, 1L))
    moved <- moved + !identical(structureOf(sampler), at)
  }
  moved
}
frozenProbs <- mixtureControl(
  proposal.probs = c(
    birth_death = 0,
    swap = 0,
    change = 0,
    perturb = 0,
    rule_gibbs = 0
  )
)@proposal.probs
# every move proposed, none at its default probability
movingProbs <- c(
  birth_death = 0.3,
  swap = 0.1,
  change = 0.2,
  perturb = 0.3,
  rule_gibbs = 0.1,
  birth = 0.4
)

# ---- a frozen mixture holds both forests -------------------------------------

held <- heteroscedastic(mixtureControl(), 2L)
set.seed(3L)
invisible(held$run(50L, 1L))
# the reader sees a variance forest that has grown and still moves
expect_true(any(varianceStructure(held)$vars >= 0L))
expect_true(sweepsMoved(held, varianceStructure) > 0L)

frozenControl <- held$control
frozenControl@proposal.probs <- frozenProbs
held$setControl(frozenControl)
surfaceAtInstall <- held$getVariance()
meanAtInstall <- meanStructure(held)
expect_identical(sweepsMoved(held, varianceStructure), 0L)
expect_identical(meanStructure(held), meanAtInstall)
# the leaf factors are still drawn
expect_false(identical(held$getVariance(), surfaceAtInstall))

# a copy re-creates the engine from the stored control: it and the sampler it
# was taken from run the same kernel
duplicate <- held$copy()
expect_identical(varianceStructure(duplicate), varianceStructure(held))
expect_identical(sweepsMoved(duplicate, varianceStructure), 0L)
expect_identical(sweepsMoved(held, varianceStructure), 0L)

# thawed by setControl, the variance forest moves again
thawedControl <- held$control
thawedControl@proposal.probs <- mixtureControl()@proposal.probs
held$setControl(thawedControl)
expect_true(sweepsMoved(held, varianceStructure) > 0L)
rm(held, duplicate)

# ---- installed before any run, it is creation under the control --------------

createdMoving <- heteroscedastic(
  mixtureControl(proposal.probs = movingProbs),
  5L
)
installedMoving <- heteroscedastic(mixtureControl(), 5L)
movingControl <- installedMoving$control
movingControl@proposal.probs <- createdMoving$control@proposal.probs
installedMoving$setControl(movingControl)
atMoving <- drawsOf(createdMoving, 6L, 20L, 20L)
expect_identical(drawsOf(installedMoving, 6L, 20L, 20L), atMoving)
expect_false(
  identical(
    drawsOf(heteroscedastic(mixtureControl(), 5L), 6L, 20L, 20L),
    atMoving
  )
)

# the width alone, under a mixture that perturbs in both forests
createdWide <- heteroscedastic(
  mixtureControl(proposal.probs = movingProbs, n.perturb.cuts = 3L),
  5L
)
installedWide <- heteroscedastic(
  mixtureControl(proposal.probs = movingProbs),
  5L
)
wideControl <- installedWide$control
wideControl@n.perturb.cuts <- 3
installedWide$setControl(wideControl)
atWide <- drawsOf(createdWide, 6L, 20L, 20L)
expect_identical(drawsOf(installedWide, 6L, 20L, 20L), atWide)
expect_false(identical(atWide, atMoving))

# ---- nothing moves for a sampler that changes neither ------------------------

# its own control handed back, and its own model: the mixture re-installed is
# the one in force on both forests
untouched <- heteroscedastic(mixtureControl(proposal.probs = movingProbs), 7L)
reinstalled <- heteroscedastic(mixtureControl(proposal.probs = movingProbs), 7L)
atUntouched <- drawsOf(untouched, 8L, 20L, 10L)
expect_identical(drawsOf(reinstalled, 8L, 20L, 10L), atUntouched)
reinstalled$setControl(reinstalled$control)
reinstalled$setModel(reinstalled$model)
expect_identical(
  drawsOf(reinstalled, 9L, 0L, 10L),
  drawsOf(untouched, 9L, 0L, 10L)
)

# ---- every chain, a variance forest's own tree prior, a single tree ----------

# the variance trees' shapes and split variables, chain by chain
chainStructures <- function(sampler) {
  sampler$storeState()
  lapply(sampler$state, function(chainState) {
    list(vars = chainState$variance.vars, sizes = chainState$variance.sizes)
  })
}
chainsMoved <- function(sampler, sweeps) {
  at <- chainStructures(sampler)
  moved <- integer(length(at))
  for (sweep in seq_len(sweeps)) {
    invisible(sampler$run(0L, 1L))
    moved <- moved + !mapply(identical, chainStructures(sampler), at)
  }
  moved
}
freeze <- function(sampler) {
  control <- sampler$control
  control@proposal.probs <- frozenProbs
  sampler$setControl(control)
  invisible(NULL)
}

# the install reaches the variance forest of every chain
set.seed(12L)
twoChains <- dbarts::dbarts(
  x,
  y,
  variance = varianceForest(n.trees = 5L),
  control = dbarts::dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE
  )
)
set.seed(13L)
invisible(twoChains$run(50L, 1L))
expect_true(all(chainsMoved(twoChains, 30L) > 0L))
freeze(twoChains)
expect_identical(chainsMoved(twoChains, 30L), c(0L, 0L))
rm(twoChains)

# and a variance forest of one tree
set.seed(14L)
oneTree <- dbarts::dbarts(
  x,
  y,
  variance = varianceForest(n.trees = 1L),
  control = mixtureControl()
)
set.seed(15L)
invisible(oneTree$run(50L, 1L))
expect_true(chainsMoved(oneTree, 100L) > 0L)
freeze(oneTree)
expect_identical(chainsMoved(oneTree, 100L), 0L)
rm(oneTree)

# a variance forest that states its own tree prior keeps it: installed before
# any run, the mixture gives creation under that prior, which is not creation
# under the mean forest's
ownPrior <- function(control, seed) {
  set.seed(seed)
  dbarts::dbarts(
    x,
    y,
    variance = varianceForest(n.trees = 5L, base = 0.7, power = 1.2),
    control = control
  )
}
createdOwn <- ownPrior(mixtureControl(proposal.probs = movingProbs), 5L)
installedOwn <- ownPrior(mixtureControl(), 5L)
ownControl <- installedOwn$control
ownControl@proposal.probs <- createdOwn$control@proposal.probs
installedOwn$setControl(ownControl)
atOwn <- drawsOf(createdOwn, 6L, 20L, 20L)
expect_identical(drawsOf(installedOwn, 6L, 20L, 20L), atOwn)
expect_false(identical(atOwn, atMoving))
