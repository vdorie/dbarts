# Tests for the per-observation ("partial") rollback mode of setPredictor:
#   sampler$setPredictor(x, column, forceUpdate = "partial")
# installs each observation's new value individually and rolls back only those whose
# value would empty a leaf, returning a per-observation logical of what was installed.
#
# A "partial" update never changes tree structure (it only re-routes observations), so a
# sampler can be returned to a known column with another partial update; the tests use this
# to restore state between scenarios rather than deep-copying (a deep copy does not preserve
# the live tree state).

# reconstruct, for a single-tree / single-predictor sampler, the (lo, hi] interval of each
# leaf from getTrees() (flattened pre-order; var == -1 marks a leaf, x > value goes right)
source(
  system.file("common", "treeLeafIntervals.R", package = "dbarts"),
  local = TRUE
)


## ---------------------------------------------------------------------------
## basic contract: a no-op proposal installs everything and changes nothing
## ---------------------------------------------------------------------------
set.seed(1L)
n <- 30L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 10L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  seed = 7L
)
sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(sampler$run(20L, 1L))

orig <- as.numeric(sampler$data@x[, 1L])
installed <- sampler$setPredictor(orig, 1L, forceUpdate = "partial")
expect_true(is.logical(installed))
expect_equal(length(installed), n)
expect_true(all(installed))
expect_equal(as.numeric(sampler$data@x[, 1L]), orig)

rm(sampler, x, y, orig, installed, n, ctrl)


## ---------------------------------------------------------------------------
## single observation: the per-observation decision (and the resulting column)
## must exactly match the existing whole-vector single-column rollback, on both
## the install and the reject branch.
## ---------------------------------------------------------------------------
set.seed(42L)
n <- 40L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 1L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  seed = 7L
)
sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(sampler$run(40L, 1L))

col <- 1L
orig <- as.numeric(sampler$data@x[, col])
leaves <- parseLeaves(sampler$getTrees())
leafOf <- leafOfObservations(leaves, orig)
sizes <- vapply(leaves, function(l) l$n, numeric(1))

srcLeaf <- which(sizes >= 2)[1]
expect_true(!is.na(srcLeaf))
members <- which(leafOf == srcLeaf)
m <- length(members)
otherObs <- which(leafOf != srcLeaf)[1] # an observation in a different, non-empty leaf

# install branch: moving one member of a multi-occupant leaf out keeps the leaf non-empty
j <- members[1]
xnew <- orig
xnew[j] <- orig[otherObs]
instP <- sampler$setPredictor(xnew, col, forceUpdate = "partial")
colP <- as.numeric(sampler$data@x[, col])
expect_true(instP[j])
expect_true(all(instP[-j]))
invisible(sampler$setPredictor(orig, col, forceUpdate = "partial")) # restore
okW <- isTRUE(sampler$setPredictor(xnew, col)) # whole-vector single-column
colW <- as.numeric(sampler$data@x[, col])
expect_true(okW)
expect_equal(colP, colW)
invisible(sampler$setPredictor(orig, col, forceUpdate = "partial")) # restore

# reject branch: reduce the leaf to a single occupant, then move that survivor out
reduced <- orig
reduced[members[1:(m - 1L)]] <- orig[otherObs]
instReduce <- sampler$setPredictor(reduced, col, forceUpdate = "partial")
expect_true(all(instReduce[members[1:(m - 1L)]]))
survivor <- members[m]
curCol <- as.numeric(sampler$data@x[, col])
xnew2 <- curCol
xnew2[survivor] <- orig[otherObs]
instP2 <- sampler$setPredictor(xnew2, col, forceUpdate = "partial")
colP2 <- as.numeric(sampler$data@x[, col])
expect_false(instP2[survivor])
expect_true(all(instP2[-survivor]))
okW2 <- isTRUE(sampler$setPredictor(xnew2, col)) # whole-vector rolls back
colW2 <- as.numeric(sampler$data@x[, col])
expect_false(okW2)
expect_equal(colP2, colW2) # both unchanged (== reduced)

rm(
  sampler,
  x,
  y,
  n,
  ctrl,
  col,
  orig,
  leaves,
  leafOf,
  sizes,
  srcLeaf,
  members,
  m,
  otherObs,
  j,
  xnew,
  instP,
  colP,
  okW,
  colW,
  reduced,
  instReduce,
  survivor,
  curCol,
  xnew2,
  instP2,
  colP2,
  okW2,
  colW2
)


## ---------------------------------------------------------------------------
## random whole-column proposals (chained, structure-preserving): every result
## is valid (no empty leaf) and data-consistent (new where installed, old where
## rolled back). Some proposals are partially rolled back, so it is not vacuous.
## ---------------------------------------------------------------------------
set.seed(11L)
n <- 60L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 30L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  seed = 3L
)
sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(sampler$run(30L, 1L))

set.seed(123L)
anyRolledBack <- FALSE
for (trial in 1:25) {
  before <- as.numeric(sampler$data@x[, 1L])
  xnew <- rnorm(n)
  inst <- sampler$setPredictor(xnew, 1L, forceUpdate = "partial")
  after <- as.numeric(sampler$data@x[, 1L])

  expect_equal(after[inst], xnew[inst])
  expect_equal(after[!inst], before[!inst])

  trees <- sampler$getTrees()
  expect_false(any(trees$var == -1 & trees$n == 0)) # no empty leaf

  if (any(!inst)) anyRolledBack <- TRUE
}
expect_true(anyRolledBack)

rm(sampler, x, y, n, ctrl, before, xnew, inst, after, trees, anyRolledBack)


## ---------------------------------------------------------------------------
## interleaving partial updates with run() (multi-chain): run() reshapes the
## trees each iteration, growing small leaves that split on the updated column
## -- the realistic outer-sampler regime. An installed observation must never
## empty a leaf, so the LIVE trees stay valid throughout. keepTrees = FALSE, so
## getTrees() reports the working trees (saved snapshots re-counted against a
## changed predictor can legitimately show empty leaves and are not the live
## state). Locks the finalize() empty-leaf invariant across every chain.
## ---------------------------------------------------------------------------
set.seed(321L)
n <- 80L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 2L,
  n.trees = 40L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  keepTrees = FALSE,
  verbose = FALSE,
  seed = 22L
)
sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(sampler$run(50L, 1L))

cur <- as.numeric(sampler$data@x[, 1L])
anyRolledBack <- FALSE
# 8 reps, not 30: n.chains = 2 makes every run() call thread-pool-bound, and at
# sd = 1.5 every one of the original 30 reps already rolled back an observation
# (checked by replay), so 8 keeps the interleave/rollback contract with margin.
for (trial in 1:8) {
  xnew <- cur + rnorm(n, sd = 1.5) # aggressive proposal exercises rollback
  inst <- sampler$setPredictor(xnew, 1L, forceUpdate = "partial")
  cur[inst] <- xnew[inst]

  trees <- sampler$getTrees() # live trees, all chains
  expect_false(any(trees$var == -1 & trees$n == 0)) # no empty leaf, any chain

  if (any(!inst)) {
    anyRolledBack <- TRUE
  }

  invisible(sampler$run(0L, 1L)) # reshape trees before the next update
}
expect_true(anyRolledBack) # the rollback path was actually exercised

rm(sampler, x, y, n, ctrl, cur, xnew, inst, trees, anyRolledBack, trial)


## ---------------------------------------------------------------------------
## coupled observations (single tree): emptying a leaf is the only constraint,
## so moving every observation out of one leaf installs all but one of them (a
## sequential sweep keeps the last occupant) regardless of the randomized scan
## order, and which one is rolled back varies across calls.
## ---------------------------------------------------------------------------
set.seed(5L)
n <- 50L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 1L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  seed = 9L
)
sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(sampler$run(40L, 1L))

col <- 1L
orig <- as.numeric(sampler$data@x[, col])
leaves <- parseLeaves(sampler$getTrees())
leafOf <- leafOfObservations(leaves, orig)

srcCandidates <- which(vapply(leaves, function(l) l$n, numeric(1)) >= 2)
expect_true(length(srcCandidates) > 0L) # the scenario must actually arise
srcLeaf <- srcCandidates[1]
members <- which(leafOf == srcLeaf)
m <- length(members)
targetObs <- which(leafOf != srcLeaf)[1]
xnew <- orig
xnew[members] <- orig[targetObs] # route every member out of the source leaf

rejected <- integer(0)
for (rep in 1:40) {
  inst <- sampler$setPredictor(xnew, col, forceUpdate = "partial")
  expect_equal(sum(inst[members]), m - 1L) # exactly one member rolled back
  expect_true(all(inst[-members])) # everyone else installed
  trees <- sampler$getTrees()
  expect_false(any(trees$var == -1 & trees$n == 0))
  rejected <- c(rejected, members[!inst[members]])
  invisible(sampler$setPredictor(orig, col, forceUpdate = "partial")) # structure-preserving restore
}
# the randomized scan order should not always reject the same member
expect_true(length(unique(rejected)) > 1L)

rm(
  sampler,
  x,
  y,
  n,
  ctrl,
  col,
  orig,
  leaves,
  leafOf,
  srcCandidates,
  srcLeaf,
  members,
  m,
  targetObs,
  xnew,
  rejected,
  inst,
  trees
)


## ---------------------------------------------------------------------------
## argument checking
## ---------------------------------------------------------------------------
set.seed(2L)
n <- 20L
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- x1 + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 5L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  seed = 1L
)
sampler <- dbarts::dbarts(
  y ~ x1 + x2,
  data.frame(x1 = x1, x2 = x2, y = y),
  control = ctrl
)
invisible(sampler$run(5L, 1L))

# partial requires a single column
expect_error(
  sampler$setPredictor(x2, forceUpdate = "partial"),
  "partial updates require a single 'column' to be specified"
)
expect_error(
  sampler$setPredictor(
    cbind(x1, x2),
    c(1L, 2L),
    forceUpdate = "partial"
  ),
  "partial updates can only be applied to a single column"
)
# partial cannot also update cut points
expect_error(
  sampler$setPredictor(
    x2,
    1L,
    forceUpdate = "partial",
    updateCutPoints = TRUE
  ),
  "partial updates cannot also update cut points"
)
# a partial update needs one value per observation, not per some other count
expect_error(
  sampler$setPredictor(x2[1:5], "x1", forceUpdate = "partial"),
  "requires one value per observation"
)
# referring to the column by name works
inst <- sampler$setPredictor(x2, "x1", forceUpdate = "partial")
expect_equal(length(inst), n)

rm(sampler, x1, x2, y, n, ctrl, inst)
rm(parseLeaves, leafOfObservations)
