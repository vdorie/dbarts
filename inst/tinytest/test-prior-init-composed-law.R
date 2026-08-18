# The no-empty-leaf conditioning the initializer applies
# (docs/design/empty-leaf-veto.md, test-prior-init-empty-leaves.R) is PER
# FOREST. A move vetoes forest f's trees against the COMPOSED precisions - the
# coupling's own per-forest weights, which carry the glue, times whatever
# weight is installed for f - so each forest's rejection draw must condition on
# exactly that vector. The chain's global weights are the wrong one by
# DEFAULT, not in a corner: a two-forest construction seeds its amplitudes at
# (a, b0, b1) = (1, 0, 1), which leaves every control row weightless in the
# treatment forest.
#
# The oracle is the sibling file's - route ONLY the rows the forest's vector
# reaches through the drawn trees (getTrees(newdata = )) and read the per-node
# counts, so a leaf no such row reaches reports n == 0 - taken through the
# internal per-forest reader, the public $getTrees having no forest axis.

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
# also make the wrong law's illegal leaves common rather than rare.
deep <- list(base = 0.95, power = 0.5)
treePrior <- dbarts::dbartsPriors$cgm(power = deep$power, base = deep$base)

makeBCF <- function(...) {
  dbarts::dbarts(
    x,
    y,
    control = control,
    tree.prior = treePrior,
    forests = list(
      dbarts::forest(),
      dbarts::forest(
        basis = ~ factor(z),
        base = deep$base,
        power = deep$power,
        ...
      )
    )
  )
}

# the internal per-forest tree reader; its forest index is 0-based (0
# prognostic, 1 treatment). rows, when given, are routed through the drawn
# trees so 'n' counts THEM per node rather than the training rows.
forestNodes <- function(sampler, forest, rows = NULL) {
  dbarts:::bartcoreGetTrees(
    list(ptr = sampler$getPointer(), x = x),
    chainNums = 1L,
    treeNums = seq_len(numTrees),
    current = TRUE,
    newdata = if (is.null(rows)) NULL else x[rows, , drop = FALSE],
    forest = forest
  )
}
leavesPerTree <- function(nodes) tapply(nodes$var == -1L, nodes$tree, sum)

# --- the per-forest law, on the DEFAULT construction ------------------------
# One pass collects both claims: no treatment-forest leaf that no treated row
# reaches, and the tree SHAPE that conditioning on 40 rows rather than 80
# produces. The shape matters because support is not the whole of a law: the
# composed conditioning removes whole trees, so it tilts the leaf-count
# distribution away from the global law's, which a single forest over the same
# design, prior and cut grid draws from.
sampler <- makeBCF()
expect_equal(as.vector(sampler$getForestAmplitudes()), c(1, 0, 1))
reference <- dbarts::dbarts(x, y, control = control, tree.prior = treePrior)

tauLeaves <- 0L
tauUnreached <- 0L
muUnreached <- 0L
tauCounts <- numeric(0)
referenceCounts <- numeric(0)
for (i in seq_len(150L)) {
  sampler$sampleTreesFromPrior()
  tau <- forestNodes(sampler, 1L, treated)
  mu <- forestNodes(sampler, 0L, treated)
  tauLeaves <- tauLeaves + sum(tau$var == -1L)
  tauUnreached <- tauUnreached + sum(tau$var == -1L & tau$n == 0L)
  muUnreached <- muUnreached + sum(mu$var == -1L & mu$n == 0L)
  tauCounts <- c(tauCounts, leavesPerTree(tau))
  reference$sampleTreesFromPrior()
  referenceCounts <- c(
    referenceCounts,
    leavesPerTree(forestNodes(reference, 0L))
  )
}
expect_true(tauLeaves > 5000L) # the trees grew; the check bites
expect_equal(tauUnreached, 0L)
# non-vacuity, and the control arm of the same law: the prognostic forest
# reaches every row (a = 1), so a leaf only control rows reach is LEGAL there
# and the same routing over the same draws finds plenty
expect_true(muUnreached > 200L)
# the shape gap, measured over 3000 trees per arm: 2.97 leaves per treatment
# tree against 3.46 for the global law, a gap of 0.48 with a Monte Carlo error
# near 0.09, and 0.023 more single-leaf trees, error near 0.009. Conditioning
# on the global vector instead closes both to noise (measured -0.05 and 0.001),
# so the bars sit between the two rather than on either.
expect_true(mean(referenceCounts) - mean(tauCounts) > 0.2)
expect_true(mean(tauCounts == 1) - mean(referenceCounts == 1) > 0.005)

# --- the empty conditioning event: bare roots, never a fault ----------------
# With no row carrying positive weight the conditional law does not exist, and
# the bare root is the only structure a later weight restore cannot strand a
# member-empty leaf in - every row sits in its one leaf. Reachable per forest,
# through a zero per-forest weight ...
zeroed <- makeBCF()
zeroed$setForestWeights(2L, rep(0, n))
zeroed$sampleTreesFromPrior()
tauNodes <- forestNodes(zeroed, 1L)
expect_equal(nrow(tauNodes), numTrees) # one node per tree
expect_true(all(tauNodes$var == -1L))
expect_true(all(tauNodes$value == 0))
# ... while the OTHER forest of the same sampler draws its own law untouched
expect_true(nrow(forestNodes(zeroed, 0L)) > 2L * numTrees)
zeroedDraws <- zeroed$run(2L, 2L)
expect_true(all(is.finite(zeroedDraws$train)))

# ... and globally, with no coupling at all: an all-zero active-row mask
# composes into the chain's own weights, a state a host whose stratum has
# emptied reaches and one the sampler accepts and runs
# (test-active-rows-pins.R). Drawing from the prior in it used to exhaust the
# rejection cap and fault.
masked <- dbarts::dbarts(x, y, control = control, tree.prior = treePrior)
masked$setActiveRows(rep(0, n))
masked$sampleTreesFromPrior()
bare <- masked$getTrees(current = TRUE)
expect_equal(nrow(bare), numTrees)
expect_true(all(bare$var == -1L))
expect_true(all(bare$value == 0))
maskedDraws <- masked$run(2L, 2L)
expect_true(all(is.finite(maskedDraws$train)))
# the guard is on the event, not on the sampler: with the mask lifted the same
# sampler draws grown trees again
masked$setActiveRows(rep(1, n))
masked$sampleTreesFromPrior()
expect_true(nrow(masked$getTrees(current = TRUE)) > 2L * numTrees)

# --- grow-from-root holds the same law -------------------------------------
# XBART-style initialization keeps both children non-empty through the scan's
# occupancy gate, which counted MEMBERS while the moves count positive weight,
# so it could land exactly the states the prior draw is now conditioned away
# from. The amplitude is held fixed so the composed vector - w * b_z^2 * s - is
# the same one across every repetition.
grown <- makeBCF(update.amplitude = FALSE)
forestWeight <- as.double(x[, 1L] <= 0.5)
grown$setForestWeights(2L, forestWeight)
reached <- treated & forestWeight > 0
expect_true(sum(reached) > 5L && sum(reached) < n %/% 3L)

grownLeaves <- 0L
grownUnreached <- 0L
for (i in seq_len(30L)) {
  grown$growFromRoot(1L)
  tau <- forestNodes(grown, 1L, reached)
  grownLeaves <- grownLeaves + sum(tau$var == -1L)
  grownUnreached <- grownUnreached + sum(tau$var == -1L & tau$n == 0L)
}
expect_true(grownLeaves > 30L * numTrees) # the trees grew; the check bites
expect_equal(grownUnreached, 0L)

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
  tauUnreached,
  muUnreached,
  tauCounts,
  referenceCounts,
  i,
  tau,
  mu,
  zeroed,
  tauNodes,
  zeroedDraws,
  masked,
  bare,
  maskedDraws,
  grown,
  forestWeight,
  reached,
  grownLeaves,
  grownUnreached
)
