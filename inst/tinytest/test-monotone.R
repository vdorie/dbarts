# Per-variable monotone (mBART) constraints: the
# `monotone` surface, its refusals, and family/prior forcing, plus a recovery
# smoke test. The exact-posterior gate lives in benchmarks/R/monotone-reference.R.

set.seed(22L)

# ---- surface: direction resolution and forcing ----

monotoneOf <- function(...) {
  sampler <- dbarts::dbarts(
    ...,
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L),
    family = gaussian(sigma = fixed(1))
  )
  probs <- sampler$control@proposal.probs
  list(
    directions = attr(sampler$model, "monotone"),
    p.birth_death = probs[["birth_death"]],
    p.swap = probs[["swap"]],
    p.change = probs[["change"]],
    p.perturb = probs[["perturb"]],
    p.rule_gibbs = probs[["rule_gibbs"]],
    k = sampler$model@leaf.hyperprior
  )
}

n <- 60L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] - x[, 3L] + rnorm(n, 0, 0.2)

# words and integers resolve to {-1, 0, +1} by column name
expect_equal(
  monotoneOf(x, y, monotone = c(a = "increasing", c = "decreasing"))$directions,
  c(1L, 0L, -1L)
)
expect_equal(
  monotoneOf(x, y, monotone = c(a = 1L, c = -1L))$directions,
  c(1L, 0L, -1L)
)
# a named 0 is unconstrained, as in the positional form
expect_equal(
  monotoneOf(x, y, monotone = c(a = 1, b = 0, c = -1))$directions,
  c(1L, 0L, -1L)
)
# an unnamed length-p vector is positional, 0 unconstrained; words and codes
# may share one, which c() makes a character vector
expect_equal(
  monotoneOf(x, y, monotone = c(1L, 0L, -1L))$directions,
  c(1L, 0L, -1L)
)
expect_equal(
  monotoneOf(x, y, monotone = c("increasing", 0, -1))$directions,
  c(1L, 0L, -1L)
)
# matching is exact and case-sensitive, as match.arg's is: the sign glyphs, a
# case variant, an abbreviation and a non-code number are refused, naming
# the vocabulary
for (bad in list("+", "-", "Increasing", "DECREASING", "inc", 2, 0.5, NA)) {
  expect_error(
    monotoneOf(x, y, monotone = list(a = bad)),
    "monotone directions must be one of \"increasing\", \"decreasing\", 1",
    fixed = TRUE
  )
}
rm(bad)

# monotone() carries the directions and the prior; the plain vector is
# shorthand for it at the default prior
spec <- dbarts::dbartsForests$monotone(c(a = "increasing"))
expect_inherits(spec, "dbartsMonotone")
expect_identical(spec$prior, "leaf")
expect_identical(
  spec,
  dbarts::dbartsForests$monotone(c(a = "increasing"), prior = "leaf")
)
expect_identical(
  monotoneOf(x, y, monotone = monotone(c(a = "increasing")))$directions,
  monotoneOf(x, y, monotone = c(a = "increasing"))$directions
)
# it has no print method of its own; the default print shows both parts
printed <- capture.output(print(spec))
expect_true(any(grepl("increasing", printed, fixed = TRUE)))
expect_true(any(grepl("leaf", printed, fixed = TRUE)))
expect_error(
  dbarts::dbartsForests$monotone(c(a = 1), prior = "tree"),
  "should be one of"
)
expect_error(dbarts::dbartsForests$monotone(), "requires 'directions'")
# the vocabulary is checked at fit time, where the columns resolve
expect_error(
  monotoneOf(x, y, monotone = monotone(c(a = "+"), prior = "leaf")),
  "must be one of"
)
# the joint prior waits for the move that samples it exactly
expect_error(
  monotoneOf(x, y, monotone = monotone(c(a = 1), prior = "joint")),
  "arrives with the corrected birth/death move"
)

# monotone() resolves by bare name inside the argument, also forwarded
# through a wrapper's dots (monotoneOf) and on dbartsSpec; a bare name the
# caller bound is that value, while a call is still the constructor
expect_equal(
  monotoneOf(x, y, monotone = monotone(c(c = "decreasing")))$directions,
  c(0L, 0L, -1L)
)
expect_equal(
  attr(
    dbarts::dbartsSpec(
      dbarts::dbartsData(x, y),
      control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L),
      monotone = monotone(c(b = "increasing"))
    )$model,
    "monotone"
  ),
  c(0L, 1L, 0L)
)
local({
  monotone <- c(a = "decreasing")
  expect_equal(monotoneOf(x, y, monotone = monotone)$directions, c(-1L, 0L, 0L))
  expect_equal(
    monotoneOf(x, y, monotone = monotone(monotone))$directions,
    c(-1L, 0L, 0L)
  )
})
# the multinomial door resolves it before refusing it
expect_error(
  dbarts::bart(
    x,
    cbind(rpois(n, 2), rpois(n, 2)),
    family = "multinomial",
    monotone = monotone(c(a = "increasing")),
    verbose = FALSE
  ),
  "does not support 'monotone'"
)

# a predictor named prior is constrained through the directions
xPrior <- x
colnames(xPrior)[2L] <- "prior"
expect_equal(
  monotoneOf(
    xPrior,
    y,
    monotone = monotone(c(prior = "increasing"), prior = "leaf")
  )$directions,
  c(0L, 1L, 0L)
)
expect_equal(
  monotoneOf(xPrior, y, monotone = c(prior = "decreasing"))$directions,
  c(0L, -1L, 0L)
)
rm(spec, printed, xPrior)

# a monotone fit is forced to birth/death-only, fixed k = 2
forced <- monotoneOf(x, y, monotone = c(a = "increasing"))
expect_equal(forced$p.birth_death, 1)
expect_equal(forced$p.swap, 0)
expect_equal(forced$p.change, 0)
expect_equal(forced$p.perturb, 0)
expect_equal(forced$p.rule_gibbs, 0)
expect_inherits(forced$k, "dbartsFixedHyperprior")
expect_equal(forced$k@k, 2)

# no constraint leaves the model unmarked and the default proposals
plain <- monotoneOf(x, y)
expect_null(plain$directions)
expect_equal(plain$p.birth_death, 0.6)

# an all-zero spec is treated as no constraint
expect_null(monotoneOf(x, y, monotone = c(a = 0L, b = 0L, c = 0L))$directions)

# ---- refusals ----

# a categorical predictor cannot carry a direction
xf <- data.frame(a = runif(n), g = factor(sample(letters[1:3], n, TRUE)))
yf <- xf$a + rnorm(n, 0, 0.2)
expect_error(
  dbarts::dbarts(yf ~ ., data = xf, monotone = c(g = "increasing")),
  "categorical"
)

# an unrecognized name is an error
expect_error(
  monotoneOf(x, y, monotone = c(nosuchcolumn = "increasing")),
  "unrecognized"
)

# an explicit k hyperprior conflicts with the fixed-k monotone rule
expect_error(
  dbarts::dbarts(
    x,
    y,
    monotone = c(a = "increasing"),
    leaf.prior = normal(chi(1.5, 2)),
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
  ),
  "monotone"
)

# an explicit non-default proposal.probs conflicts with birth/death-only
expect_error(
  dbarts::dbarts(
    x,
    y,
    monotone = c(a = "increasing"),
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0.1,
      change = 0.3,
      birth = 0.5
    ),
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
  ),
  "proposal.probs"
)

# but a caller spelling the documented default and omitting the move that
# ships at zero is not a non-default vector: the refusal fills the name it
# compares before reading it, so this is rewritten rather than refused
expect_equal(
  monotoneOf(
    x,
    y,
    monotone = c(a = "increasing"),
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.4,
      birth = 0.5
    )
  )$p.birth_death,
  1
)

# an explicit proposal.probs = NULL is treated as absent: it succeeds
# under monotone and is bitwise identical to the defaulted call
ctrlD4 <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 5L,
  seed = 55L
)
defaultedD4 <- dbarts::dbarts(
  x,
  y,
  monotone = c(a = "increasing"),
  control = ctrlD4
)
defaultedD4$sampleTreesFromPrior()
samplesDefaultedD4 <- defaultedD4$run(0L, 5L)

explicitNullD4 <- dbarts::dbarts(
  x,
  y,
  monotone = c(a = "increasing"),
  proposal.probs = NULL,
  control = ctrlD4
)
explicitNullD4$sampleTreesFromPrior()
samplesExplicitNullD4 <- explicitNullD4$run(0L, 5L)

expect_identical(samplesDefaultedD4$train, samplesExplicitNullD4$train)
expect_identical(samplesDefaultedD4$sigma, samplesExplicitNullD4$sigma)

rm(
  ctrlD4,
  defaultedD4,
  samplesDefaultedD4,
  explicitNullD4,
  samplesExplicitNullD4
)

rm(monotoneOf, forced, plain, x, y, xf, yf)

# ---- recovery: a monotone truth is recovered with tighter intervals ----

set.seed(101L)
nRec <- 200L
xRec <- matrix(sort(runif(nRec)), nRec, 1L, dimnames = list(NULL, "x1"))
truth <- 1.5 * xRec[, 1L]
yRec <- truth + rnorm(nRec, 0, 0.3)

fitArgs <- list(
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 500L,
  n.burn = 300L,
  keepTrees = FALSE,
  verbose = FALSE
)
fitMono <- do.call(
  dbarts::bart,
  c(list(xRec, yRec, monotone = c(x1 = "increasing")), fitArgs)
)
fitFree <- do.call(dbarts::bart, c(list(xRec, yRec), fitArgs))

# the constrained posterior-mean fit is monotone in x (x is sorted)
expect_true(all(diff(fitMono$yhat.train.mean) > -1e-8))

# it recovers the truth
expect_true(
  sqrt(mean((fitMono$yhat.train.mean - truth)^2)) / sd(truth) < 0.25
)

# and its pointwise posterior intervals are, on average, tighter than the
# unconstrained fit's (the constraint pools information across the ordering)
sdMono <- mean(apply(fitMono$yhat.train, 2L, sd))
sdFree <- mean(apply(fitFree$yhat.train, 2L, sd))
expect_true(sdMono < sdFree)

# ---- a non-monotone truth flattens under the constraint, does not crash ----

set.seed(202L)
yBump <- sin(2 * pi * xRec[, 1L]) + rnorm(nRec, 0, 0.3)
fitBump <- do.call(
  dbarts::bart,
  c(list(xRec, yBump, monotone = c(x1 = "increasing")), fitArgs)
)
fitBumpFree <- do.call(dbarts::bart, c(list(xRec, yBump), fitArgs))
# the fit stays monotone (the true rise-then-fall is flattened, not recovered)
expect_true(all(diff(fitBump$yhat.train.mean) > -1e-8))
expect_true(all(is.finite(fitBump$yhat.train.mean)))
# the flattening bias: the unconstrained fit recovers the rise-and-fall (large
# range); the monotone fit cannot chase the descent, so its range is smaller
expect_true(
  diff(range(fitBump$yhat.train.mean)) <
    diff(range(fitBumpFree$yhat.train.mean))
)

rm(
  nRec,
  xRec,
  truth,
  yRec,
  fitArgs,
  fitMono,
  fitFree,
  sdMono,
  sdFree,
  yBump,
  fitBump,
  fitBumpFree
)

# ---- the prior draw is constrained too ----

# samplePriorPredictive installs its draws as the sampler's live leaf values,
# so under a constraint they must come from the truncated prior; an
# unconstrained draw would both misstate the prior and strand the sampler in an
# infeasible state.

set.seed(303L)
nPri <- 80L
xPri <- matrix(sort(runif(nPri)), nPri, 1L, dimnames = list(NULL, "x1"))
samplerPri <- dbarts::dbarts(
  xPri,
  xPri[, 1L] + rnorm(nPri, 0, 0.3),
  monotone = c(x1 = "increasing"),
  control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 25L)
)

priorEv <- samplePriorPredictive(samplerPri, n.samples = 100L, type = "ev")
expect_equal(dim(priorEv), c(100L, nPri))
expect_true(all(is.finite(priorEv)))

# x is sorted, so each individual draw - not merely their average - is
# non-decreasing across the constrained column
expect_true(all(apply(priorEv, 1L, function(f) all(diff(f) > -1e-8))))

# the prior draw is not the constant zero forest it would collapse to if the
# rejection simply gave up
expect_true(mean(apply(priorEv, 1L, function(f) diff(range(f)))) > 0.1)

# and a sampler left in that state runs on
samplerPri$sampleTreesFromPrior()
samplerPri$sampleLeafParametersFromPrior()
samplesPri <- samplerPri$run(10L, 10L)
expect_true(all(is.finite(samplesPri$train)))

rm(nPri, xPri, samplerPri, priorEv, samplesPri)

# ---- reachability: every state the sampler holds lies in the cone ----

# x1 carries the constraint and x2 is free; predictions along x1 at fixed x2
# read the fitted function, including leaves a test row reaches but no
# training row does
set.seed(404L)
nReach <- 200L
xReach <- matrix(
  runif(nReach * 2L),
  nReach,
  2L,
  dimnames = list(NULL, c("x1", "x2"))
)
xReach[1L, ] <- 0
xReach[2L, ] <- 1
# decreasing along the constrained axis, so the unconstrained fit and the
# constrained likelihood both pull against the cone
yReach <- -2 * xReach[, 1L] + sin(4 * xReach[, 2L]) + rnorm(nReach, 0, 0.2)
controlReach <- function(n.trees = 20L, ...) {
  dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    n.samples = 1L,
    updateState = TRUE,
    ...
  )
}
gridReach <- as.matrix(expand.grid(
  x1 = seq(0, 1, length.out = 101L),
  x2 = c(0.1, 0.5, 0.9)
))
# the largest fall along x1 at any fixed x2; 0 for a monotone fit
maxDrop <- function(sampler) {
  fits <- matrix(sampler$predict(gridReach), 101L)
  -min(apply(fits, 2L, diff))
}
monoReach <- dbarts::dbarts(
  xReach,
  yReach,
  monotone = c(x1 = "increasing"),
  control = controlReach()
)
freeReach <- dbarts::dbarts(xReach, yReach, control = controlReach())
invisible(freeReach$run(200L, 1L))
invisible(monoReach$run(50L, 1L))
birthDeath <- c(
  birth_death = 1,
  swap = 0,
  change = 0,
  perturb = 0,
  rule_gibbs = 0,
  birth = 0.5
)

# setControl mirrors creation: a defaulted mixture is rewritten to
# birth/death-only, and a mixture proposing other moves is refused with the
# control left as it was
monoReach$setControl(controlReach(printEvery = 50L))
expect_equal(monoReach$control@proposal.probs, birthDeath)
expect_error(
  monoReach$setControl(controlReach(
    proposal.probs = c(birth_death = 0.5, change = 0.5)
  )),
  "birth/death-only"
)
expect_equal(monoReach$control@proposal.probs, birthDeath)

# the bridge refuses the mixture on its own, keyed on the engine's leaf kind:
# a control slipped past the R check meets it at the prior install
changeControl <- monoReach$control
changeControl@proposal.probs[] <- c(0.5, 0, 0.5, 0, 0, 0.5)
storedControl <- monoReach$control
monoReach$control <- changeControl
expect_error(
  monoReach$setModel(monoReach$model),
  "proposes only birth and death"
)
monoReach$control <- storedControl
rm(changeControl, storedControl)

# a warm start from an unconstrained donor installs its trees, reseeding
# every tree outside the cone to all-zero, so the fit is monotone at once
expect_true(maxDrop(freeReach) > 0.1)
monoReach$installTrees(freeReach)
expect_true(maxDrop(monoReach) <= 1e-8)
invisible(monoReach$run(0L, 1L))
expect_true(maxDrop(monoReach) <= 1e-8)

# setState of the donor's state is refused whole, and the sampler keeps its
# state and its fit
stateBefore <- monoReach$state
fitBefore <- monoReach$predict(gridReach)
expect_error(
  monoReach$setState(freeReach$state),
  "leaf values violate this sampler's monotone constraint"
)
expect_identical(monoReach$state, stateBefore)
expect_identical(monoReach$predict(gridReach), fitBefore)
rm(stateBefore, fitBefore)

# grow-from-root regrows every tree over recycled node slots; the regrown
# leaves draw from the all-zero seed, and one sweep later the fit is monotone
monoReach$growFromRoot(2L)
expect_true(maxDrop(monoReach) <= 1e-8)
invisible(monoReach$run(0L, 1L))
expect_true(maxDrop(monoReach) <= 1e-8)

# a long run against the constraint stays in the cone: a copy (setState
# through the stored state) and a state round trip both install
invisible(monoReach$run(2000L, 1L))
expect_true(maxDrop(monoReach) <= 1e-8)
monoCopy <- monoReach$copy()
expect_identical(monoCopy$predict(gridReach), monoReach$predict(gridReach))
monoReach$setState(monoReach$state)
expect_true(maxDrop(monoCopy) <= 1e-8)
rm(monoCopy)

# a hand-built tree over x1 (increasing) and x2 (free): x1 splits at 0.5 and
# each side splits x2 at 0.3, leaves J, K (x1 low) and S1, S2 (x1 high), with
# J < S1 and K < S2 the only relations. A forced update emptying S1 collapses
# x1's high side into one leaf holding S2's value, which now borders J from
# above and sits below it; the merged tree is reseeded, so the fit stays
# monotone and the state reinstalls
handReach <- dbarts::dbarts(
  xReach,
  yReach,
  monotone = c(x1 = "increasing"),
  control = controlReach(n.trees = 1L)
)
handState <- handReach$state
cuts <- attr(handState, "cutPoints")
cut1 <- cuts[[1L]][which.max(cuts[[1L]] >= 0.5)]
cut2 <- cuts[[2L]][which.max(cuts[[2L]] >= 0.3)]
handTree <- handState[[1L]]$forests[[1L]]
handTree$tree.vars <- c(1L, 2L, -1L, -1L, 2L, -1L, -1L)
handTree$tree.values <- writeBin(
  c(cut1, cut2, 0, -0.1, cut2, 0.02, -0.09),
  raw()
)
handTree$tree.sizes <- 7L
handTree$tree.flags <- as.raw(c(2L, 2L, 0L, 0L, 2L, 0L, 0L))
handState[[1L]]$forests[[1L]] <- handTree
handReach$setState(handState)
xEmpty <- xReach
movedRows <- xEmpty[, 1L] > cut1 & xEmpty[, 2L] <= cut2
xEmpty[movedRows, 2L] <- (cut2 + 1) / 2
expect_true(any(movedRows))
handReach$setPredictor(xEmpty, forceUpdate = TRUE)
expect_true(maxDrop(handReach) <= 1e-8)
handReach$storeState()
expect_silent(handReach$setState(handReach$state))
expect_equal(handReach$state[[1L]]$forests[[1L]]$tree.sizes, 5L)
rm(handReach, handState, cuts, cut1, cut2, handTree, xEmpty, movedRows)

# a forced predictor update that flattens x2 empties one side of every x2
# split; the collapse merges leaves across the free axis, relating each merged
# leaf to neighbors it did not border, and a merged tree outside the cone is
# reseeded
xFlat <- xReach
xFlat[-(1:2), 2L] <- 0.5
monoReach$setPredictor(xFlat, forceUpdate = TRUE)
expect_true(maxDrop(monoReach) <= 1e-8)
monoReach$storeState()
monoReach$setState(monoReach$state)
invisible(monoReach$run(20L, 1L))
expect_true(maxDrop(monoReach) <= 1e-8)
rm(xFlat)

# empty leaves carry drawn values in the cone: a 20-tree fit warm-started
# from a monotone donor over [0, 1] into rows that leave x1 in (0.6, 1) empty
# (the endpoints keep the grids equal, so the trees install verbatim) strands
# leaves where the fit sits above 0. The structure is then frozen (the
# all-zero mixture, which a monotone sampler allows), so the empty leaves are
# redrawn every sweep instead of being merged away by the first death. Every
# tree is checked on its own, read through getTrees and evaluated at the
# midpoint of every cell its splits cut out, so no other tree's rise can mask
# a fall and no cell is missed
donorReach <- dbarts::dbarts(
  xReach,
  -yReach,
  monotone = c(x1 = "increasing"),
  control = controlReach()
)
invisible(donorReach$run(300L, 1L))
xGap <- xReach
xGap[-(1:2), 1L] <- runif(nReach - 2L, 0, 0.6)
gapReach <- dbarts::dbarts(
  xGap,
  -yReach,
  monotone = c(x1 = "increasing"),
  control = controlReach()
)
gapReach$installTrees(donorReach)
gapReach$setControl(controlReach(
  proposal.probs = c(birth_death = 0, change = 0)
))
# one tree's leaf values at x, its nodes in pre-order (var -1 marks a leaf, a
# row goes left when x[, var] <= value)
treeAt <- function(nodes, x) {
  position <- 0L
  descend <- function(rows) {
    position <<- position + 1L
    var <- nodes$var[position]
    value <- nodes$value[position]
    if (var < 0L) {
      return(rep(value, length(rows)))
    }
    goesLeft <- x[rows, var] <= value
    out <- numeric(length(rows))
    out[goesLeft] <- descend(rows[goesLeft])
    out[!goesLeft] <- descend(rows[!goesLeft])
    out
  }
  descend(seq_len(nrow(x)))
}
# the largest fall along x1 in any single tree, over every cell of that tree
worstTreeDrop <- function(sampler) {
  trees <- sampler$getTrees()
  drops <- vapply(
    split(trees, trees$tree),
    function(nodes) {
      mids <- lapply(1:2, function(j) {
        edges <- sort(unique(c(0, 1, nodes$value[nodes$var == j])))
        (edges[-1L] + edges[-length(edges)]) / 2
      })
      cells <- as.matrix(expand.grid(mids[[1L]], mids[[2L]]))
      fits <- matrix(treeAt(nodes, cells), length(mids[[1L]]))
      if (nrow(fits) < 2L) 0 else -min(apply(fits, 2L, diff))
    },
    0.0
  )
  max(drops)
}
# the sum of the trees read this way is the sampler's own fit up to its
# response scale, so the reader routes as the engine does
treesGap <- gapReach$getTrees()
totalGap <- Reduce(
  `+`,
  lapply(split(treesGap, treesGap$tree), treeAt, x = gridReach)
)
expect_equal(cor(totalGap, gapReach$predict(gridReach)), 1)
expect_true(any(treesGap$var < 0L & treesGap$n == 0L))
dropsGap <- vapply(
  seq_len(20L),
  function(i) {
    invisible(gapReach$run(0L, 1L))
    worstTreeDrop(gapReach)
  },
  0.0
)
expect_true(all(dropsGap <= 1e-12))
treesGap <- gapReach$getTrees()
expect_true(any(treesGap$var < 0L & treesGap$n == 0L))

rm(
  nReach,
  xReach,
  yReach,
  controlReach,
  gridReach,
  maxDrop,
  monoReach,
  freeReach,
  birthDeath,
  donorReach,
  xGap,
  gapReach,
  treeAt,
  worstTreeDrop,
  treesGap,
  totalGap,
  dropsGap
)
