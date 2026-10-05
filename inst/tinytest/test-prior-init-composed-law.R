# The no-empty-leaf conditioning the initializer applies
# (test-prior-init-empty-leaves.R) is ONE law for every forest. A move refuses
# a tree of forest f that leaves a leaf no row reaches, and reads neither the
# coupling's per-forest precisions - which carry the glue - nor a weight
# installed for f, so each forest's rejection draw conditions on membership
# alone. That matters by DEFAULT, not in a corner: a two-forest construction
# seeds its amplitudes at (a, b0, b1) = (1, 0, 1), which leaves every control
# row weightless in the treatment forest, and a treatment tree may still hold a
# leaf of control rows only.
#
# The oracle is the sibling file's - the drawn trees report the rows that reach
# each node, and routing ONLY the rows a forest's precisions reach through them
# (getTrees(newdata = )) the rows its likelihood sees - taken through the
# public $getTrees' own forest argument.

set.seed(20260818L)
n <- 80L
x <- matrix(runif(n * 2L), n, 2L)
z <- rep(0:1, length.out = n)
treated <- z == 1L
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

numTrees <- 20L
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = numTrees,
  updateState = FALSE,
  seed = 11L
)
# a deep growth prior, stated for BOTH forests so the treatment forest does not
# fall back to the shallower per-forest default: the conditioning is what these
# arms measure, and a difference in the prior would stand in for it. Deep trees
# also make leaves of few rows common rather than rare.
deep <- list(base = 0.95, power = 0.5)
treePrior <- dbarts::dbartsPriors$cgm(power = deep$power, base = deep$base)

makeBCF <- function(...) {
  dbarts::dbarts(
    x,
    y,
    control = control,
    tree.prior = treePrior,
    forests = list(
      forest(),
      forest(
        basis = ~ factor(z),
        base = deep$base,
        power = deep$power,
        ...
      )
    )
  )
}

# the per-forest tree reader; its forest index is 1-based (1 prognostic, 2
# treatment). rows, when given, are routed through the drawn trees so 'n'
# counts THEM per node rather than the training rows.
forestNodes <- function(sampler, forest, rows = NULL) {
  sampler$getTrees(
    forest = forest,
    treeNums = seq_len(numTrees),
    chainNums = 1L,
    current = TRUE,
    newdata = if (is.null(rows)) NULL else x[rows, , drop = FALSE]
  )
}
leavesPerTree <- function(nodes) tapply(nodes$var == -1L, nodes$tree, sum)

# --- one law for both forests, on the DEFAULT construction ------------------
# One pass collects three claims: every leaf of either forest holds a row; a
# treatment-forest leaf that no treated row reaches is legal and common; and
# the treatment trees have the SHAPE the same prior draws over all 80 rows in a
# single forest of the same design, prior and cut grid, not the smaller one
# that conditioning on the 40 treated rows would give.
sampler <- makeBCF()
expect_equal(as.vector(sampler$getForestAmplitudes()), c(1, 0, 1))
reference <- dbarts::dbarts(x, y, control = control, tree.prior = treePrior)

tauLeaves <- 0L
tauEmpty <- 0L
muEmpty <- 0L
tauUnreached <- 0L
tauCounts <- numeric(0)
referenceCounts <- numeric(0)
for (i in seq_len(150L)) {
  sampler$sampleTreesFromPrior()
  tauMembers <- forestNodes(sampler, 2L)
  muMembers <- forestNodes(sampler, 1L)
  tau <- forestNodes(sampler, 2L, treated)
  tauLeaves <- tauLeaves + sum(tau$var == -1L)
  tauEmpty <- tauEmpty + sum(tauMembers$var == -1L & tauMembers$n == 0L)
  muEmpty <- muEmpty + sum(muMembers$var == -1L & muMembers$n == 0L)
  tauUnreached <- tauUnreached + sum(tau$var == -1L & tau$n == 0L)
  tauCounts <- c(tauCounts, leavesPerTree(tau))
  reference$sampleTreesFromPrior()
  referenceCounts <- c(
    referenceCounts,
    leavesPerTree(forestNodes(reference, 1L))
  )
}
expect_true(tauLeaves > 5000L) # the trees grew; the check bites
expect_equal(tauEmpty, 0L)
expect_equal(muEmpty, 0L)
expect_true(tauUnreached > 200L)
# the shape, over 3000 trees per arm: conditioning the treatment forest on its
# 40 treated rows alone would leave it 0.48 leaves per tree short of the
# single forest's 3.46, against a Monte Carlo error near 0.05; the two are
# 0.01 apart
expect_true(abs(mean(referenceCounts) - mean(tauCounts)) < 0.2)

# --- no row weighted at all: the same law, never a fault ---------------------
# A forest whose every row carries zero weight draws its trees from the same
# conditioned prior, since the conditioning does not read the weights.
# Reachable per forest, through a zero per-forest weight ...
zeroed <- makeBCF()
zeroed$setForestWeights(2L, rep(0, n))
zeroed$sampleTreesFromPrior()
tauNodes <- forestNodes(zeroed, 2L)
expect_true(nrow(tauNodes) > 2L * numTrees)
expect_true(all(tauNodes$n[tauNodes$var == -1L] > 0L))
expect_true(all(tauNodes$value[tauNodes$var == -1L] == 0))
expect_true(nrow(forestNodes(zeroed, 1L)) > 2L * numTrees)
zeroedDraws <- zeroed$run(2L, 2L)
expect_true(all(is.finite(zeroedDraws$train)))

# ... and globally, with no coupling at all: an all-zero active-row mask is a
# state a host whose stratum has emptied reaches and one the sampler accepts
# and runs (test-active-rows-pins.R)
masked <- dbarts::dbarts(x, y, control = control, tree.prior = treePrior)
masked$setActiveRows(rep(0, n))
masked$sampleTreesFromPrior()
maskedNodes <- masked$getTrees(current = TRUE)
expect_true(nrow(maskedNodes) > 2L * numTrees)
expect_true(all(maskedNodes$n[maskedNodes$var == -1L] > 0L))
maskedDraws <- masked$run(2L, 2L)
expect_true(all(is.finite(maskedDraws$train)))
# with the mask lifted the same sampler draws from the same law
masked$setActiveRows(rep(1, n))
masked$sampleTreesFromPrior()
expect_true(nrow(masked$getTrees(current = TRUE)) > 2L * numTrees)

# --- grow-from-root holds the same law -------------------------------------
# XBART-style initialization keeps both children non-empty through the scan's
# occupancy gate, which counts MEMBERS as the moves do: every grown leaf holds
# a row, and one may hold only rows the forest's precisions do not reach. The
# amplitude is held fixed so the composed vector - w * b_z^2 * s - is the same
# one across every repetition.
grown <- makeBCF(update.amplitude = FALSE)
forestWeight <- as.double(x[, 1L] <= 0.5)
grown$setForestWeights(2L, forestWeight)
reached <- treated & forestWeight > 0
expect_true(sum(reached) > 5L && sum(reached) < n %/% 3L)

grownLeaves <- 0L
grownEmpty <- 0L
grownUnreached <- 0L
for (i in seq_len(30L)) {
  grown$growFromRoot(1L)
  tauMembers <- forestNodes(grown, 2L)
  tau <- forestNodes(grown, 2L, reached)
  grownLeaves <- grownLeaves + sum(tau$var == -1L)
  grownEmpty <- grownEmpty + sum(tauMembers$var == -1L & tauMembers$n == 0L)
  grownUnreached <- grownUnreached + sum(tau$var == -1L & tau$n == 0L)
}
expect_true(grownLeaves > 30L * numTrees) # the trees grew; the check bites
expect_equal(grownEmpty, 0L)
expect_true(grownUnreached > 0L)

rm(
  n,
  x,
  z,
  treated,
  y,
  numTrees,
  control,
  deep,
  treePrior,
  makeBCF,
  forestNodes,
  leavesPerTree,
  sampler,
  reference,
  tauLeaves,
  tauEmpty,
  muEmpty,
  tauUnreached,
  tauCounts,
  referenceCounts,
  i,
  tau,
  tauMembers,
  muMembers,
  zeroed,
  tauNodes,
  zeroedDraws,
  masked,
  maskedNodes,
  maskedDraws,
  grown,
  forestWeight,
  reached,
  grownLeaves,
  grownEmpty,
  grownUnreached
)
