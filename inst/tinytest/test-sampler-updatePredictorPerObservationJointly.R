# Tests for the joint per-observation rollback mode:
#   updatePredictorPerObservationJointly(list(samplerA, samplerB, ...), x, column)
# A single sequential sweep installs each observation in EVERY sampler at once, and only if its
# new value keeps every leaf non-empty in every tree of every chain of every sampler. Returns the
# per-observation install mask (identical for all samplers); afterwards every sampler holds the
# same values in the shared column.
#
# Like the single-fit "partial" mode, a joint update never changes tree structure (it only
# re-routes observations), so the samplers can be returned to a known shared column with another
# joint update -- the tests use that to restore state between scenarios.

# reconstruct, for a single-tree / single-predictor sampler, the (lo, hi] interval of each leaf
# from getTrees() (flattened pre-order; var == -1 marks a leaf, x > value goes right)
parseLeaves <- function(rows) {
  i <- 0L
  leaves <- list()
  recurse <- function(lo, hi) {
    i <<- i + 1L
    row <- rows[i, ]
    if (row$var == -1) {
      leaves[[length(leaves) + 1L]] <<- list(lo = lo, hi = hi, n = row$n)
      return(invisible())
    }
    v <- row$value
    recurse(lo, v) # left child:  x <= value
    recurse(v, hi) # right child: x >  value
  }
  recurse(-Inf, Inf)
  leaves
}
leafOfObservations <- function(leaves, values) {
  vapply(
    values,
    function(xj) {
      which(vapply(leaves, function(l) l$lo < xj & xj <= l$hi, logical(1)))
    },
    integer(1)
  )
}

# single-predictor sampler on a shared x with its own response (so trees differ across samplers)
makeSampler <- function(y, x, seed, ntree = 1L) {
  ctrl <- dbarts::dbartsControl(
    n.chains = 1L,
    n.trees = ntree,
    n.samples = 1L,
    n.burn = 0L,
    updateState = FALSE,
    verbose = FALSE,
    rngSeed = seed
  )
  s <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
  invisible(s$run(40L, 1L))
  s
}
noEmptyLeaf <- function(s) {
  trees <- s$getTrees()
  !any(trees$var == -1 & trees$n == 0)
}


## ---------------------------------------------------------------------------
## no-op: installing the current (shared) column installs everything and leaves
## both columns unchanged and identical.
## ---------------------------------------------------------------------------
set.seed(101L)
n <- 50L
x <- rnorm(n)
A <- makeSampler(x + rnorm(n), x, 11L)
B <- makeSampler(-x + rnorm(n), x, 22L)

orig <- as.numeric(A$data@x[, 1L])
expect_equal(orig, as.numeric(B$data@x[, 1L])) # the column is genuinely shared

inst <- updatePredictorPerObservationJointly(list(A, B), orig, "x")
expect_true(is.logical(inst))
expect_equal(length(inst), n)
expect_true(all(inst))
expect_equal(as.numeric(A$data@x[, 1L]), orig)
expect_equal(as.numeric(B$data@x[, 1L]), orig)

rm(A, B, x, n, orig, inst)


## ---------------------------------------------------------------------------
## single-observation change: with only one observation moving there is no
## coupling, so the scan order is irrelevant and the joint decision is exactly
## the AND of the two single-fit "partial" decisions:
##   instJoint[j] == instA[j] && instB[j],   instJoint[-j] all TRUE.
## This is the defining semantics; we check it for every observation, moving each
## into a different leaf of A so the reject branch (and the role of B) is exercised.
## ---------------------------------------------------------------------------
set.seed(202L)
n <- 60L
x <- rnorm(n)
A <- makeSampler(x + rnorm(n), x, 7L)
B <- makeSampler(sin(3 * x) + rnorm(n), x, 99L)
orig <- as.numeric(A$data@x[, 1L])

leavesA <- parseLeaves(A$getTrees())
leafA <- leafOfObservations(leavesA, orig)
for (j in seq_len(n)) {
  otherLeafObs <- which(leafA != leafA[j])
  if (length(otherLeafObs) == 0L) {
    next
  } # A is a single leaf: nothing to move into
  xnew <- orig
  xnew[j] <- orig[otherLeafObs[1L]]

  # both samplers start each scenario at the shared 'orig'
  instA <- A$setPredictor(xnew, "x", forceUpdate = "partial")
  invisible(A$setPredictor(orig, "x", forceUpdate = "partial")) # restore A
  instB <- B$setPredictor(xnew, "x", forceUpdate = "partial")
  invisible(B$setPredictor(orig, "x", forceUpdate = "partial")) # restore B

  instJ <- updatePredictorPerObservationJointly(list(A, B), xnew, "x")

  expect_equal(instJ[j], instA[j] && instB[j]) # the joint decision is the AND
  expect_true(all(instJ[-j])) # untouched observations always install

  cA <- as.numeric(A$data@x[, 1L])
  cB <- as.numeric(B$data@x[, 1L])
  expect_equal(cA, cB) # identical after the joint call
  expect_equal(cA[j], if (instJ[j]) xnew[j] else orig[j])
  expect_true(noEmptyLeaf(A))
  expect_true(noEmptyLeaf(B))

  invisible(updatePredictorPerObservationJointly(list(A, B), orig, "x")) # structure-preserving restore in both
}

rm(
  A,
  B,
  x,
  n,
  orig,
  leavesA,
  leafA,
  j,
  otherLeafObs,
  xnew,
  instA,
  instB,
  instJ,
  cA,
  cB
)


## ---------------------------------------------------------------------------
## coupled observations: routing EVERY member of one (single-tree) sampler's leaf
## out at once cannot empty that leaf, so the joint sweep must keep at least one
## member back -- a guaranteed rejection that also keeps both columns identical and
## leaves no empty leaf in either sampler. (The single-tree A is the binding fit; B
## only adds constraints, so the rejection count is never below one.)
## ---------------------------------------------------------------------------
set.seed(303L)
n <- 50L
x <- rnorm(n)
A <- makeSampler(x + rnorm(n), x, 13L) # single tree -> clean leaf-interval oracle
B <- makeSampler(exp(x / 2) + rnorm(n), x, 31L, ntree = 15L)
orig <- as.numeric(A$data@x[, 1L])

leavesA <- parseLeaves(A$getTrees())
leafA <- leafOfObservations(leavesA, orig)
sizes <- vapply(leavesA, function(l) l$n, numeric(1))
srcLeaf <- which(sizes >= 2L)[1L]
expect_true(!is.na(srcLeaf)) # the scenario must arise
members <- which(leafA == srcLeaf)
m <- length(members)
targetObs <- which(leafA != srcLeaf)[1L]
xnew <- orig
xnew[members] <- orig[targetObs] # route every member of A's leaf out

inst <- updatePredictorPerObservationJointly(list(A, B), xnew, "x")
expect_true(sum(!inst[members]) >= 1L) # at least one member is rolled back
expect_true(all(inst[-members])) # everyone else installs
cA <- as.numeric(A$data@x[, 1L])
cB <- as.numeric(B$data@x[, 1L])
expect_equal(cA, cB)
expect_true(noEmptyLeaf(A))
expect_true(noEmptyLeaf(B))

rm(
  A,
  B,
  x,
  n,
  orig,
  leavesA,
  leafA,
  sizes,
  srcLeaf,
  members,
  m,
  targetObs,
  xnew,
  inst,
  cA,
  cB
)


## ---------------------------------------------------------------------------
## random whole-column proposals (multi-tree, chained, structure-preserving):
## after every joint call the two samplers' shared columns are identical, are
## data-consistent (new where installed, old where rolled back), and neither has
## an empty leaf. Some proposals roll observations back, so it is not vacuous.
## ---------------------------------------------------------------------------
set.seed(303L)
n <- 80L
x <- rnorm(n)
A <- makeSampler(x + rnorm(n), x, 5L, ntree = 20L)
B <- makeSampler(cos(2 * x) + rnorm(n), x, 6L, ntree = 20L)

set.seed(404L)
anyRolledBack <- FALSE
for (trial in 1:30) {
  before <- as.numeric(A$data@x[, 1L])
  expect_equal(before, as.numeric(B$data@x[, 1L]))
  xnew <- rnorm(n)
  inst <- updatePredictorPerObservationJointly(list(A, B), xnew, "x")
  cA <- as.numeric(A$data@x[, 1L])
  cB <- as.numeric(B$data@x[, 1L])

  expect_equal(cA, cB) # the invariant the old resync loop had to force
  expect_equal(cA[inst], xnew[inst])
  expect_equal(cA[!inst], before[!inst])
  expect_true(noEmptyLeaf(A))
  expect_true(noEmptyLeaf(B))

  if (any(!inst)) anyRolledBack <- TRUE
}
expect_true(anyRolledBack)

rm(A, B, x, n, before, xnew, inst, cA, cB, anyRolledBack, trial)


## ---------------------------------------------------------------------------
## three-way joint with the shared column at different positions in each sampler
## (resolved by name), plus argument checking.
## ---------------------------------------------------------------------------
set.seed(505L)
n <- 40L
x <- rnorm(n)
w <- rnorm(n)
v <- rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 10L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE,
  rngSeed = 3L
)
A <- dbarts::dbarts(
  y ~ theta + w,
  data.frame(theta = x, w = w, y = x + rnorm(n)),
  control = ctrl
)
B <- dbarts::dbarts(
  y ~ v + theta,
  data.frame(v = v, theta = x, y = -x + rnorm(n)),
  control = ctrl
)
invisible(A$run(20L, 1L))
invisible(B$run(20L, 1L))

# theta sits in column 1 of A but column 2 of B -- resolution is by name
expect_equal(match("theta", colnames(A$data@x)), 1L)
expect_equal(match("theta", colnames(B$data@x)), 2L)

xnew <- rnorm(n)
inst <- updatePredictorPerObservationJointly(list(A, B), xnew, "theta")
expect_equal(length(inst), n)
expect_equal(as.numeric(A$data@x[, "theta"]), as.numeric(B$data@x[, "theta"]))
expect_equal(as.numeric(A$data@x[inst, "theta"]), xnew[inst])
expect_true(noEmptyLeaf(A))
expect_true(noEmptyLeaf(B))

# a single sampler may be passed directly (not wrapped in a list)
inst2 <- updatePredictorPerObservationJointly(
  A,
  as.numeric(A$data@x[, "theta"]),
  "theta"
)
expect_true(all(inst2))

# argument checking
expect_error(updatePredictorPerObservationJointly(list(A, B), xnew)) # column missing
expect_error(updatePredictorPerObservationJointly(
  list(A, B),
  xnew,
  c("theta", "w")
)) # not a single column
expect_error(updatePredictorPerObservationJointly(list(A, 1), xnew, "theta")) # non-sampler in samplers
expect_error(updatePredictorPerObservationJointly(list(), xnew, "theta")) # empty sampler list
C <- dbarts::dbarts(
  y ~ q,
  data.frame(q = rnorm(n), y = rnorm(n)),
  control = ctrl
)
invisible(C$run(5L, 1L))
expect_error(updatePredictorPerObservationJointly(list(A, C), xnew, "theta")) # column absent in a sampler
D <- dbarts::dbarts(
  y ~ theta,
  data.frame(theta = rnorm(n + 5L), y = rnorm(n + 5L)),
  control = ctrl
)
invisible(D$run(5L, 1L))
expect_error(updatePredictorPerObservationJointly(list(A, D), xnew, "theta")) # observation count mismatch

rm(A, B, C, D, x, w, v, n, ctrl, xnew, inst, inst2)


## ---------------------------------------------------------------------------
## interleaving joint updates with run() (multi-chain): run() reshapes each
## sampler's trees between updates -- the realistic outer-sampler regime. An
## installed observation must keep every leaf non-empty in BOTH samplers across
## every chain, so the live trees stay valid and the two shared columns stay
## identical throughout. keepTrees = FALSE, so getTrees() reports the working
## trees. Locks the finalize() empty-leaf invariant on the joint path.
## ---------------------------------------------------------------------------
set.seed(654L)
nr <- 80L
xr <- rnorm(nr)
rctrl <- dbarts::dbartsControl(
  n.chains = 2L,
  n.trees = 40L,
  n.samples = 1L,
  n.burn = 0L,
  updateState = FALSE,
  keepTrees = FALSE,
  verbose = FALSE,
  rngSeed = 44L
)
RA <- dbarts::dbarts(
  y ~ theta,
  data.frame(theta = xr, y = xr + rnorm(nr)),
  control = rctrl
)
RB <- dbarts::dbarts(
  y ~ theta,
  data.frame(theta = xr, y = -xr + rnorm(nr)),
  control = rctrl
)
invisible(RA$run(50L, 1L))
invisible(RB$run(50L, 1L))

cur <- as.numeric(RA$data@x[, "theta"])
anyRolledBack <- FALSE
# 8 reps, not 30: n.chains = 2 makes every run() call thread-pool-bound, and at
# sd = 1.5 every one of the original 30 reps already rolled back an observation
# (checked by replay), so 8 keeps the interleave/rollback contract with margin.
for (trial in 1:8) {
  prop <- cur + rnorm(nr, sd = 1.5) # aggressive proposal exercises rollback
  mask <- updatePredictorPerObservationJointly(list(RA, RB), prop, "theta")
  cur[mask] <- prop[mask]

  expect_true(noEmptyLeaf(RA)) # live trees valid, all chains
  expect_true(noEmptyLeaf(RB))
  expect_equal(
    as.numeric(RA$data@x[, "theta"]),
    as.numeric(RB$data@x[, "theta"])
  ) # shared column stays in sync

  if (any(!mask)) {
    anyRolledBack <- TRUE
  }

  invisible(RA$run(0L, 1L))
  invisible(RB$run(0L, 1L))
}
expect_true(anyRolledBack) # the rollback path was actually exercised

rm(RA, RB, nr, xr, rctrl, cur, prop, mask, anyRolledBack, trial)

rm(parseLeaves, leafOfObservations, makeSampler, noEmptyLeaf)
