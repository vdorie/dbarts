# A state goes onto the sampler's cut points and does not replace them. A
# live setState writes neither the points nor their weights: a state whose
# split thresholds are all points of the sampler's grid installs as stored,
# TRUE where its trees fit; one with a threshold the grid does not hold is
# declined, FALSE and the sampler as it was, and with forceUpdate = TRUE has
# that split moved to the nearest point its tree still allows, or merged with
# what is beneath it into one leaf at the mean of the stored leaves weighted
# by the rows each held. That is what installing the state on its own points
# and then calling setCutPoints(splits = "value") gives. A copy and a reload,
# which create a new engine, still take the points and weights their stored
# state records. To restore a state's cut points too the caller sets them
# first: setCutPoints(attr(st, "cutPoints")) and then setState(st).

visibly <- function(value) list(value = value, visible = TRUE)
forcedValue <- list(value = NULL, visible = FALSE)
# the state the engine holds now, the sampler's state field left alone
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
cutPointsOf <- function(sampler) attr(stateOf(sampler), "cutPoints")
massOf <- function(sampler) attr(stateOf(sampler), "cutMass")
treeFields <- c("tree.vars", "tree.values", "tree.sizes", "tree.flags")
# the live mean trees of every chain, as stored
liveOf <- function(state) {
  lapply(state, function(chain) lapply(chain$forests, `[`, treeFields))
}
reloadOf <- function(sampler) {
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  reloaded <- readRDS(path)
  unlink(path)
  reloaded
}
# A state with tree 1 of a chain's forest replaced by a hand-built one in
# pre-order: vars holds a column, counted from one, for a split and -1 for a
# leaf, values a threshold or a leaf value, and flags 2 for a split on a
# threshold, plus 1 where it sends a missing value right.
withTree <- function(
  state,
  vars,
  values,
  chain = 1L,
  flags = ifelse(vars > 0, 2L, 0L)
) {
  forest <- state[[chain]]$forests[[1L]]
  first <- seq_len(forest$tree.sizes[1L])
  forest$tree.vars <- c(as.integer(vars), forest$tree.vars[-first])
  forest$tree.values <- c(
    writeBin(as.double(values), raw()),
    forest$tree.values[-seq_len(8L * length(first))]
  )
  forest$tree.flags <- c(as.raw(flags), forest$tree.flags[-first])
  forest$tree.sizes[1L] <- length(vars)
  state[[chain]]$forests[[1L]] <- forest
  state
}
firstTree <- function(state, chain = 1L) {
  forest <- state[[chain]]$forests[[1L]]
  size <- forest$tree.sizes[1L]
  list(
    vars = forest$tree.vars[seq_len(size)],
    values = readBin(forest$tree.values, "double", size),
    flags = as.integer(forest$tree.flags[seq_len(size)])
  )
}
# every threshold the live mean trees of a state hold on a column
thresholdsOn <- function(state, column) {
  unlist(lapply(state, function(chain) {
    forest <- chain$forests[[1L]]
    values <- readBin(forest$tree.values, "double", length(forest$tree.vars))
    values[forest$tree.vars == column]
  }))
}
expectError <- function(expr, info) {
  expect_error(
    expr,
    pattern = "state is not consistent with this sampler",
    fixed = TRUE,
    info = info
  )
}

set.seed(11L, sample.kind = "Rejection")
n <- 200L
x <- cbind(a = runif(n), b = runif(n))
y <- sin(4 * x[, 1L]) + x[, 2L] + rnorm(n, 0, 0.2)
controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 4L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
# a sampler past a short burn-in, or one whose trees are all single leaves;
# built twice the same way, two are the same sampler
make <- function(..., run = TRUE, data = x, response = y) {
  sampler <- dbarts::dbarts(data, response, control = controlWith(...))
  if (run) {
    invisible(sampler$run(20L, 1L))
  }
  sampler
}
grid <- cutPointsOf(make(run = FALSE))
fine <- grid[[1L]]
expect_identical(lengths(grid), c(100L, 100L))
three <- c(0.25, 0.5, 0.75)
# three of the points the samplers are created with, and every one of them
# with the points halfway between
coarse <- fine[c(25L, 50L, 75L)]
finer <- sort(c(fine, (fine[-1L] + fine[-100L]) / 2))
expect_identical(length(unique(finer)), 199L)

# ---- live grid kept -------------------------------------------------------

# three points set, an earlier state forced in: the sampler holds three
# points, and the weights it held, as before
sampler <- make()
stored <- stateOf(sampler)
sampler$setCutPoints(three, 1L)
before <- stateOf(sampler)
expect_identical(lengths(attr(before, "cutPoints")), c(3L, 100L))
expect_identical(
  withVisible(sampler$setState(stored, forceUpdate = TRUE)),
  forcedValue
)
expect_identical(cutPointsOf(sampler), attr(before, "cutPoints"))
expect_null(massOf(sampler))
expect_true(all(thresholdsOn(stateOf(sampler), 1L) %in% three))
expect_false(all(thresholdsOn(stored, 1L) %in% three))

# ---- same grid exact ------------------------------------------------------

# store, five sweeps, restore: TRUE, the state the engine then holds is the
# stored one, and the next draws are those of a twin that did the same
sampler <- make()
twin <- make()
stored <- stateOf(sampler)
expect_identical(stateOf(twin), stored)
invisible(sampler$run(0L, 5L))
invisible(twin$run(0L, 5L))
expect_false(identical(liveOf(stateOf(sampler)), liveOf(stored)))
expect_identical(withVisible(sampler$setState(stored)), visibly(TRUE))
expect_identical(sampler$state, stored)
expect_identical(stateOf(sampler), stored)
expect_true(twin$setState(stored))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))

# ---- finer grid is clean --------------------------------------------------

# 199 points holding the stored 100: TRUE, the 199 kept, every split as
# stored; forced, the same trees
sampler <- make()
stored <- stateOf(sampler)
invisible(sampler$run(0L, 5L))
sampler$setCutPoints(finer, 1L, splits = "value")
expect_identical(withVisible(sampler$setState(stored)), visibly(TRUE))
expect_identical(cutPointsOf(sampler), list(finer, grid[[2L]]))
expect_identical(liveOf(stateOf(sampler)), liveOf(stored))
invisible(sampler$run(0L, 2L))
expect_identical(
  withVisible(sampler$setState(stored, forceUpdate = TRUE)),
  forcedValue
)
expect_identical(cutPointsOf(sampler), list(finer, grid[[2L]]))
expect_identical(liveOf(stateOf(sampler)), liveOf(stored))

# ---- on the grid, a leaf empty --------------------------------------------

# the finer grid, and the column changed so that a stored leaf holds no row:
# FALSE, as on the stored grid, and forced the leaf is merged on either
inner <- x[, 1L] > min(x[, 1L]) & x[, 1L] < max(x[, 1L])
collapsed <- ifelse(inner, 0.5, x[, 1L])
onFiner <- make()
onStored <- make()
onFiner$setCutPoints(finer, 1L, splits = "value")
for (sampler in list(onFiner, onStored)) {
  sampler$setPredictor(collapsed, 1L, forceUpdate = TRUE)
  held <- stateOf(sampler)
  expect_identical(withVisible(sampler$setState(stored)), visibly(FALSE))
  expect_identical(stateOf(sampler), held)
  expect_identical(
    withVisible(sampler$setState(stored, forceUpdate = TRUE)),
    forcedValue
  )
  expect_identical(cutPointsOf(sampler), attr(held, "cutPoints"))
}
expect_identical(liveOf(stateOf(onFiner)), liveOf(stateOf(onStored)))
expect_false(identical(liveOf(stateOf(onFiner)), liveOf(stored)))

# ---- off the grid declines ------------------------------------------------

# FALSE, and the sampler is the twin that never called: its grid, weights
# and trees, the state its engine stores, its state field and its next draws
sampler <- make()
twin <- make()
stored <- stateOf(sampler)
for (each in list(sampler, twin)) {
  invisible(each$run(0L, 3L))
  each$setCutPoints(three, 1L)
  each$storeState()
}
field <- sampler$state
expect_identical(withVisible(sampler$setState(stored)), visibly(FALSE))
expect_identical(sampler$state, field)
expect_identical(stateOf(sampler), stateOf(twin))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))

# ---- nearest inside the interval ------------------------------------------

# hand-built trees on a recorded grid that holds their thresholds, forced
# onto {0.25, 0.5, 0.75}. By position the three roots would all go to 0.75.
recorded <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.37, 0.375, 0.5, 0.6)
movedTree <- function(vars, values, ...) {
  sampler <- make(n.trees = 5L, run = FALSE)
  state <- withTree(stateOf(sampler), vars, values, ...)
  attr(state, "cutPoints")[[1L]] <- recorded
  sampler$setCutPoints(three, 1L)
  expect_identical(withVisible(sampler$setState(state)), visibly(FALSE))
  expect_identical(
    withVisible(sampler$setState(state, forceUpdate = TRUE)),
    forcedValue
  )
  expect_identical(cutPointsOf(sampler)[[1L]], three)
  firstTree(stateOf(sampler))
}
# 0.37 lies nearer 0.25 than 0.5
tree <- movedTree(c(1, -1, -1), c(0.37, 0.1, -0.2))
expect_identical(tree$vars, c(1L, -1L, -1L))
expect_identical(tree$values, c(0.25, 0.1, -0.2))
# 0.375 lies midway and goes to the upper point
tree <- movedTree(c(1, -1, -1), c(0.375, 0.1, -0.2))
expect_identical(tree$values, c(0.5, 0.1, -0.2))
# a root at 0.6 goes to 0.5, and its left child, stored on 0.5, to the
# nearest point below it and not to the root's own
tree <- movedTree(c(1, 1, -1, -1, -1), c(0.6, 0.5, 0.1, -0.2, 0.05))
expect_identical(tree$vars, c(1L, 1L, -1L, -1L, -1L))
expect_identical(tree$values, c(0.5, 0.25, 0.1, -0.2, 0.05))

# ---- merged leaf: left, mirror, two deep ----------------------------------

# a split with no point left inside its interval becomes, with everything
# beneath it, one leaf at the mean of the stored leaves weighted by the rows
# each held under the stored thresholds
mergedTree <- function(onto, vars, values, data = x, state = NULL) {
  sampler <- make(n.trees = 5L, run = FALSE, data = data)
  if (is.null(state)) {
    state <- stateOf(sampler)
  }
  state <- withTree(state, vars, values)
  sampler$setCutPoints(onto, 1L)
  expect_identical(withVisible(sampler$setState(state)), visibly(FALSE))
  expect_identical(
    withVisible(sampler$setState(state, forceUpdate = TRUE)),
    forcedValue
  )
  firstTree(stateOf(sampler))
}
# the mean of leaf values weighted by the rows of a column in each interval
# (lower, upper]
heldMean <- function(column, lower, upper, values) {
  rows <- vapply(
    seq_along(values),
    function(k) sum(column > lower[k] & column <= upper[k]),
    0L
  )
  expect_true(all(rows > 0L) && length(unique(rows)) > 1L)
  sum(rows * values) / sum(rows)
}
# left: the root on the 50th point and its left child on the 30th, onto the
# 50th and 75th
tree <- mergedTree(
  fine[c(50L, 75L)],
  c(1, 1, -1, -1, -1),
  c(fine[50L], fine[30L], 0.10, -0.20, 0.05)
)
expected <- heldMean(
  x[, 1L],
  c(-Inf, fine[30L]),
  fine[c(30L, 50L)],
  c(0.10, -0.20)
)
expect_identical(tree$vars, c(1L, -1L, -1L))
expect_identical(tree$values[c(1L, 3L)], c(fine[50L], 0.05))
expect_equal(tree$values[2L], expected, tolerance = 1e-12)
# neither the plain mean nor the leaf the rows land in once the child moved
expect_true(abs(expected + 0.05) > 0.01 && abs(expected - 0.10) > 0.01)
# mirror: the right child on the 60th point, onto the 30th and 50th
tree <- mergedTree(
  fine[c(30L, 50L)],
  c(1, -1, 1, -1, -1),
  c(fine[50L], 0.05, fine[60L], 0.10, -0.20)
)
expected <- heldMean(
  x[, 1L],
  fine[c(50L, 60L)],
  c(fine[60L], Inf),
  c(0.10, -0.20)
)
expect_identical(tree$vars, c(1L, -1L, -1L))
expect_identical(tree$values[1L:2L], c(fine[50L], 0.05))
expect_equal(tree$values[3L], expected, tolerance = 1e-12)
expect_true(abs(expected + 0.05) > 0.01 && abs(expected + 0.20) > 0.01)
# two deep: the left child and its left child both without a point
tree <- mergedTree(
  fine[c(50L, 75L)],
  c(1, 1, 1, -1, -1, -1, -1),
  c(fine[50L], fine[30L], fine[10L], 0.3, 0.1, -0.2, 0.05)
)
expected <- heldMean(
  x[, 1L],
  c(-Inf, fine[c(10L, 30L)]),
  fine[c(10L, 30L, 50L)],
  c(0.3, 0.1, -0.2)
)
expect_identical(tree$vars, c(1L, -1L, -1L))
expect_identical(tree$values[c(1L, 3L)], c(fine[50L], 0.05))
expect_equal(tree$values[2L], expected, tolerance = 1e-12)

# ---- another sampler's data -----------------------------------------------

# a state from a sampler on other rows: the stored leaves are weighed by the
# rows of the sampler the state goes into
xOther <- cbind(a = runif(n)^2, b = runif(n))
donor <- make(n.trees = 5L, run = FALSE, data = xOther)
other <- cutPointsOf(donor)[[1L]]
expect_false(identical(other, fine))
tree <- mergedTree(
  other[c(50L, 75L)],
  c(1, 1, -1, -1, -1),
  c(other[50L], other[30L], 0.10, -0.20, 0.05),
  state = stateOf(donor)
)
expected <- heldMean(
  x[, 1L],
  c(-Inf, other[30L]),
  other[c(30L, 50L)],
  c(0.10, -0.20)
)
donorMean <- heldMean(
  xOther[, 1L],
  c(-Inf, other[30L]),
  other[c(30L, 50L)],
  c(0.10, -0.20)
)
expect_true(abs(expected - donorMean) > 0.01)
expect_identical(tree$vars, c(1L, -1L, -1L))
expect_equal(tree$values[2L], expected, tolerance = 1e-12)

# ---- equals the two calls -------------------------------------------------

# a fitted state forced across grids leaves the sampler a twin is left in by
# the state installed on its own grid and setCutPoints(splits = "value")
# then called: every tree, leaf value and generator, and the next draws
donor <- make()
stored <- stateOf(donor)
sampler <- make()
invisible(sampler$run(0L, 3L))
sampler$setCutPoints(three, 1L)
expect_identical(
  withVisible(sampler$setState(stored, forceUpdate = TRUE)),
  forcedValue
)
twin <- make()
expect_true(twin$setState(stored))
twin$setCutPoints(three, 1L, splits = "value")
moved <- stateOf(sampler)
expect_identical(moved, stateOf(twin))
expect_identical(lengths(attr(moved, "cutPoints")), c(3L, 100L))
expect_true(
  length(thresholdsOn(moved, 1L)) < length(thresholdsOn(stored, 1L))
)
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))

# ---- field after crossing -------------------------------------------------

# after an install onto other cut points or weights than the state records
# the state field holds what the engine holds, so a copy and a reload have
# the sampler's grid, weights and trees
expectFieldFollows <- function(sampler, given, info) {
  held <- stateOf(sampler)
  expect_identical(sampler$state, held, info = info)
  expect_false(identical(sampler$state, given), info = info)
  for (door in list(sampler$copy(), reloadOf(sampler))) {
    remade <- stateOf(door)
    expect_identical(
      attr(remade, "cutPoints"),
      attr(held, "cutPoints"),
      info = info
    )
    expect_identical(attr(remade, "cutMass"), attr(held, "cutMass"), info)
    expect_identical(liveOf(remade), liveOf(held), info = info)
  }
}
sampler <- make()
sampler$setCutPoints(three, 1L)
expect_null(sampler$setState(stored, forceUpdate = TRUE))
expectFieldFollows(sampler, stored, "forced across grids")
sampler <- make()
sampler$setCutPoints(finer, 1L, splits = "value")
expect_true(sampler$setState(stored))
expectFieldFollows(sampler, stored, "clean on a finer grid")
# a forced repair on the state's own grid leaves the state as given
sampler <- make()
sampler$setPredictor(collapsed, 1L, forceUpdate = TRUE)
expect_null(sampler$setState(stored, forceUpdate = TRUE))
expect_identical(sampler$state, stored)
expect_false(identical(liveOf(stateOf(sampler)), liveOf(stored)))

# ---- reload, first call ---------------------------------------------------

# a sampler saved on three points takes a state stored on 100 as its first
# call after a reload exactly as after any other call: its engine is first
# rebuilt from its own saved state
sampler <- make()
invisible(sampler$run(0L, 2L))
sampler$setCutPoints(three, 1L)
sampler$storeState()
expect_false(identical(liveOf(sampler$state), liveOf(stored)))
first <- reloadOf(sampler)
later <- reloadOf(sampler)
invisible(later$getTrees())
expect_identical(withVisible(first$setState(stored)), visibly(FALSE))
expect_identical(withVisible(later$setState(stored)), visibly(FALSE))
expect_identical(stateOf(first), stateOf(later))
expect_identical(lengths(cutPointsOf(first)), c(3L, 100L))
expect_identical(first$state, sampler$state)
first <- reloadOf(sampler)
later <- reloadOf(sampler)
invisible(later$getTrees())
expect_null(first$setState(stored, forceUpdate = TRUE))
expect_null(later$setState(stored, forceUpdate = TRUE))
expect_identical(stateOf(first), stateOf(later))
expect_identical(lengths(cutPointsOf(first)), c(3L, 100L))
expect_identical(first$state, later$state)
# with no saved state of its own, or one that is refused, the new engine
# holds no grid of the sampler's and takes the handed state's
stateless <- reloadOf(sampler)
stateless$state <- NULL
expect_identical(withVisible(stateless$setState(stored)), visibly(TRUE))
expect_identical(cutPointsOf(stateless), grid)
expect_identical(stateless$state, stored)
refusing <- reloadOf(sampler)
refusing$state <- stateOf(make(n.trees = 7L))
expect_identical(withVisible(refusing$setState(stored)), visibly(TRUE))
expect_identical(cutPointsOf(refusing), grid)
expect_identical(liveOf(stateOf(refusing)), liveOf(stored))

# ---- weights -------------------------------------------------------------

# a column of 30 distinct values and one of 7 take a point in each gap,
# weighed by the gap's width
xWeighted <- cbind(
  v = rep_len(0:29, n)^2 / 100,
  s = as.double((seq_len(n) %% 7L)^2)
)
yWeighted <- sqrt(xWeighted[, 1L]) +
  xWeighted[, 2L] / 10 +
  rnorm(n, 0, 0.2)
weighted <- function(...) make(..., data = xWeighted, response = yWeighted)
sampler <- weighted()
stored <- stateOf(sampler)
points <- attr(stored, "cutPoints")
mass <- attr(stored, "cutMass")
expect_identical(lengths(points), c(29L, 6L))
expect_identical(lengths(mass), c(30L, 7L))

# the weights are the sampler's: on equal points a state recording none
# goes into a weighted sampler, and a weighted one into an unweighted, TRUE
# and forced NULL, the weights as they were and the trees as stored
bare <- stored
attr(bare, "cutMass") <- NULL
unweighted <- weighted()
for (j in 1L:2L) {
  unweighted$setCutPoints(points[[j]] + 0.125, j)
  unweighted$setCutPoints(points[[j]], j)
}
expect_identical(cutPointsOf(unweighted), points)
expect_null(massOf(unweighted))
cases <- list(
  list(sampler = sampler, state = bare, mass = mass, info = "bare state"),
  list(sampler = unweighted, state = stored, mass = NULL, info = "weighted")
)
for (case in cases) {
  into <- case$sampler
  invisible(into$run(0L, 2L))
  expect_identical(withVisible(into$setState(case$state)), visibly(TRUE))
  expect_identical(massOf(into), case$mass, info = case$info)
  expect_identical(liveOf(stateOf(into)), liveOf(stored), info = case$info)
  expect_identical(into$state, stateOf(into), info = case$info)
  invisible(into$run(0L, 2L))
  expect_identical(
    withVisible(into$setState(case$state, forceUpdate = TRUE)),
    forcedValue
  )
  expect_identical(massOf(into), case$mass, info = case$info)
  expect_identical(cutPointsOf(into), points, info = case$info)
  expect_identical(liveOf(stateOf(into)), liveOf(stored), info = case$info)
}

# by hand, half 1: three points set on the first column, then the stored
# points of the whole set. The changed column has its points and no weights,
# the unchanged one keeps its weights, the state goes in as stored and the
# weights stay as they are
sampler <- weighted()
sampler$setCutPoints(c(1, 3, 6), 1L)
expect_identical(massOf(sampler), list(NULL, mass[[2L]]))
sampler$setCutPoints(attr(stored, "cutPoints"))
expect_identical(cutPointsOf(sampler), points)
expect_identical(massOf(sampler), list(NULL, mass[[2L]]))
expect_identical(withVisible(sampler$setState(stored)), visibly(TRUE))
expect_identical(massOf(sampler), list(NULL, mass[[2L]]))
expect_identical(liveOf(stateOf(sampler)), liveOf(stored))
expectFieldFollows(sampler, stored, "other weights by hand")

# by hand, half 2: the column's own values under updateCutPoints = "value"
# give points and weights back, and the chain is the stored one
sampler <- weighted()
twin <- weighted()
sampler$setCutPoints(c(1, 3, 6), 1L)
expect_true(sampler$setPredictor(
  xWeighted[, 1L],
  1L,
  updateCutPoints = "value"
))
expect_identical(cutPointsOf(sampler), points)
expect_identical(massOf(sampler), mass)
expect_identical(withVisible(sampler$setState(stored)), visibly(TRUE))
expect_identical(sampler$state, stored)
expect_true(twin$setState(stored))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))

# ---- re-creation adopts ---------------------------------------------------

# stored after three points: a copy and a reload hold the three points and
# the weights of the other column
sampler <- weighted()
sampler$setCutPoints(c(1, 3, 6), 1L)
sampler$storeState()
for (door in list(sampler$copy(), reloadOf(sampler))) {
  expect_identical(cutPointsOf(door), list(c(1, 3, 6), points[[2L]]))
  expect_identical(massOf(door), list(NULL, mass[[2L]]))
  expect_identical(liveOf(stateOf(door)), liveOf(sampler$state))
}

# ---- by hand --------------------------------------------------------------

# after three points, with the state held in a variable, the whole stored
# set and then the state give the stored grid identically, TRUE, and the
# next draws of a twin that only restored
frame <- data.frame(
  a = x[, 1L],
  b = x[, 2L],
  f = factor(rep_len(c("p", "q", "r", "s"), n)),
  o = ordered(rep_len(c(1L, 3L, 2L, 4L, 4L), n))
)
caseWeights <- runif(n, 0.5, 2)
few <- cbind(a = as.double(seq_len(n) %% 6L), b = x[, 2L])
fixtures <- list(
  numeric = function() dbarts::dbarts(x, y, control = controlWith()),
  weights = function() {
    dbarts::dbarts(x, y, weights = caseWeights, control = controlWith())
  },
  unordered = function() {
    dbarts::dbarts(y ~ a + b + f, frame, control = controlWith())
  },
  ordered = function() {
    dbarts::dbarts(y ~ a + b + o, frame, control = controlWith())
  },
  few = function() dbarts::dbarts(few, y, control = controlWith())
)
for (name in names(fixtures)) {
  sampler <- fixtures[[name]]()
  twin <- fixtures[[name]]()
  invisible(sampler$run(20L, 1L))
  invisible(twin$run(20L, 1L))
  st <- stateOf(sampler)
  set <- if (name == "few") c(0.5, 2.5, 4.5) else three
  sampler$setCutPoints(set, 1L)
  expect_identical(cutPointsOf(sampler)[[1L]], set, info = name)
  sampler$setCutPoints(attr(st, "cutPoints"))
  expect_identical(cutPointsOf(sampler), attr(st, "cutPoints"), info = name)
  expect_null(massOf(sampler), info = name)
  expect_null(attr(st, "cutMass"), info = name)
  expect_identical(withVisible(sampler$setState(st)), visibly(TRUE), name)
  expect_identical(sampler$state, st, info = name)
  expect_true(twin$setState(st), info = name)
  expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L), info = name)
}
expect_identical(lengths(attr(st, "cutPoints")), c(5L, 100L))

# ---- refresh away and back ------------------------------------------------

# a column replaced with its grid refreshed and replaced back the same way,
# then the restore: the grid is the stored one, TRUE, and the draws are
# those of a twin that only restored
sampler <- make()
twin <- make()
stored <- stateOf(sampler)
away <- runif(n, 1, 3)
for (column in list(away, x[, 1L])) {
  sampler$setPredictor(
    column,
    1L,
    forceUpdate = TRUE,
    updateCutPoints = "position"
  )
}
expect_identical(cutPointsOf(sampler), grid)
expect_identical(withVisible(sampler$setState(stored)), visibly(TRUE))
expect_true(twin$setState(stored))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))
# put back with no refresh after one, the column's grid is the other
# values' and the state is declined
sampler <- make()
sampler$setPredictor(
  away,
  1L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
)
sampler$setPredictor(x[, 1L], 1L, forceUpdate = TRUE)
expect_false(identical(cutPointsOf(sampler), grid))
expect_identical(withVisible(sampler$setState(stored)), visibly(FALSE))

# ---- kept draws stay ------------------------------------------------------

# kept draws are copied and never moved: forced across grids, predict gives
# what the state's sampler gave and the kept thresholds are the stored ones
savedFields <- c("saved.vars", "saved.values", "saved.sizes", "saved.flags")
keeping <- function(seed) {
  sampler <- make(keepTrees = TRUE, seed = seed, run = FALSE)
  invisible(sampler$run(20L, 4L))
  sampler
}
donor <- keeping(7L)
stored <- stateOf(donor)
sampler <- keeping(8L)
sampler$setCutPoints(three, 1L)
expect_false(identical(sampler$predict(x), donor$predict(x)))
expect_null(sampler$setState(stored, forceUpdate = TRUE))
expect_identical(sampler$predict(x), donor$predict(x))
expect_identical(
  stateOf(sampler)[[1L]]$forests[[1L]][savedFields],
  stored[[1L]]$forests[[1L]][savedFields]
)
expect_identical(lengths(cutPointsOf(sampler)), c(3L, 100L))
expect_true(all(thresholdsOn(stateOf(sampler), 1L) %in% three))

# ---- variance forest ------------------------------------------------------

# a variance split off the grid, alone, declines the state; forced, it lies
# on the grid
varianceFields <- paste0("variance.", c("vars", "values", "sizes", "flags"))
varianceThresholds <- function(state, column) {
  chain <- state[[1L]]
  values <- readBin(
    chain$variance.values,
    "double",
    length(chain$variance.vars)
  )
  values[chain$variance.vars == column]
}
sampler <- dbarts::dbarts(x, y, variance = TRUE, control = controlWith())
invisible(sampler$run(20L, 1L))
sampler$setCutPoints(coarse, 1L, splits = "value")
onCoarse <- stateOf(sampler)
# the sampler's own state recorded on the 100 points, which hold its three
attr(onCoarse, "cutPoints")[[1L]] <- fine
expect_identical(withVisible(sampler$setState(onCoarse)), visibly(TRUE))
offGrid <- onCoarse
size <- offGrid[[1L]]$variance.sizes[1L]
offGrid[[1L]]$variance.vars <- c(
  c(1L, -1L, -1L),
  offGrid[[1L]]$variance.vars[-seq_len(size)]
)
offGrid[[1L]]$variance.values <- c(
  writeBin(c(fine[30L], 1.25, 0.8), raw()),
  offGrid[[1L]]$variance.values[-seq_len(8L * size)]
)
offGrid[[1L]]$variance.flags <- c(
  as.raw(c(2L, 0L, 0L)),
  offGrid[[1L]]$variance.flags[-seq_len(size)]
)
offGrid[[1L]]$variance.sizes[1L] <- 3L
expect_identical(liveOf(offGrid), liveOf(onCoarse))
held <- stateOf(sampler)
expect_identical(withVisible(sampler$setState(offGrid)), visibly(FALSE))
expect_identical(stateOf(sampler), held)
expect_null(sampler$setState(offGrid, forceUpdate = TRUE))
moved <- stateOf(sampler)
expect_identical(attr(moved, "cutPoints")[[1L]], coarse)
expect_true(all(varianceThresholds(moved, 1L) %in% coarse))
expect_identical(
  readBin(moved[[1L]]$variance.values, "double", 3L),
  c(coarse[1L], 1.25, 0.8)
)
expect_identical(liveOf(moved), liveOf(onCoarse))

# ---- second chain ---------------------------------------------------------

# an off-grid split in the second chain alone: FALSE, both chains as they
# were; forced, the first chain's trees are the stored ones
sampler <- make(n.chains = 2L)
sampler$setCutPoints(coarse, 1L, splits = "value")
onCoarse <- stateOf(sampler)
attr(onCoarse, "cutPoints")[[1L]] <- fine
expect_identical(withVisible(sampler$setState(onCoarse)), visibly(TRUE))
offGrid <- withTree(
  onCoarse,
  c(1, -1, -1),
  c(fine[30L], 0.1, -0.2),
  chain = 2L
)
expect_identical(liveOf(offGrid)[[1L]], liveOf(onCoarse)[[1L]])
invisible(sampler$run(0L, 2L))
held <- stateOf(sampler)
expect_identical(withVisible(sampler$setState(offGrid)), visibly(FALSE))
expect_identical(stateOf(sampler), held)
expect_null(sampler$setState(offGrid, forceUpdate = TRUE))
moved <- stateOf(sampler)
expect_identical(liveOf(moved)[[1L]], liveOf(onCoarse)[[1L]])
expect_identical(
  firstTree(moved, chain = 2L)$values,
  c(coarse[1L], 0.1, -0.2)
)

# ---- neither grid ---------------------------------------------------------

# a threshold that is a point of neither grid is an error, whatever grid
# the state records, and the sampler is as it was
sampler <- make()
sampler$setCutPoints(coarse, 1L, splits = "value")
held <- stateOf(sampler)
stray <- withTree(held, c(1, -1, -1), c(0.123456, 0.1, -0.2))
onOther <- stray
attr(onOther, "cutPoints")[[1L]] <- fine
# and one the recorded grid lacks though the sampler's holds it
unrecorded <- withTree(held, c(1, -1, -1), c(coarse[2L], 0.1, -0.2))
attr(unrecorded, "cutPoints")[[1L]] <- fine[-50L]
for (state in list(stray, onOther, unrecorded)) {
  expectError(sampler$setState(state), "neither grid")
  expectError(sampler$setState(state, forceUpdate = TRUE), "neither grid")
}
expect_identical(stateOf(sampler), held)

# ---- kinds ----------------------------------------------------------------

# a state records no column kind. An unordered factor's entry is empty, and
# an ordered factor's grid is its level count's, so a state recording
# anything else on such a column is another column's: an error in both
# forms, the sampler as it was
levelFrame <- function(kind) {
  codes <- rep_len(c(1L, 3L, 2L, 4L, 4L), n)
  if (kind == "five") {
    codes <- rep_len(c(1L, 3L, 2L, 5L, 4L), n)
  }
  column <- switch(
    kind,
    numeric = x[, 1L],
    whole = as.double(codes) - 1,
    unordered = factor(codes),
    ordered(codes)
  )
  data.frame(a = column, b = x[, 2L])
}
kindOf <- function(kind) {
  sampler <- dbarts::dbarts(
    y ~ a + b,
    levelFrame(kind),
    control = controlWith()
  )
  invisible(sampler$run(20L, 1L))
  sampler
}
states <- lapply(
  c(numeric = "numeric", ordered = "ordered", five = "five"),
  function(kind) stateOf(kindOf(kind))
)
expect_identical(attr(states$ordered, "cutPoints")[[1L]], c(0.5, 1.5, 2.5))
expect_identical(
  attr(states$five, "cutPoints")[[1L]],
  c(0.5, 1.5, 2.5, 3.5)
)
refusals <- list(
  list(into = "unordered", state = states$numeric),
  list(into = "ordered", state = states$numeric),
  list(into = "ordered", state = states$five)
)
for (case in refusals) {
  sampler <- kindOf(case$into)
  held <- stateOf(sampler)
  expectError(sampler$setState(case$state), case$into)
  expectError(sampler$setState(case$state, forceUpdate = TRUE), case$into)
  expect_identical(stateOf(sampler), held, info = case$into)
}
# another grid on an ordered factor is refused even where its thresholds are
# all points of the sampler's: the five-level grid holds the four-level one
sampler <- kindOf("ordered")
held <- stateOf(sampler)
relabelled <- held
attr(relabelled, "cutPoints")[[1L]] <- c(0.5, 1.5, 2.5, 3.5)
expectError(sampler$setState(relabelled), "ordered, other levels")
expectError(
  sampler$setState(relabelled, forceUpdate = TRUE),
  "ordered, other levels"
)
expect_identical(stateOf(sampler), held)
expect_true(sampler$setState(held))
# the other way round cannot be told: an ordered factor's entry is what a
# numeric column of whole numbers records, so its state follows the numeric
# rule, FALSE where a threshold is off the grid and TRUE where none is
expect_true(length(thresholdsOn(states$ordered, 1L)) > 0L)
sampler <- kindOf("numeric")
held <- stateOf(sampler)
expect_identical(withVisible(sampler$setState(states$ordered)), visibly(FALSE))
expect_identical(stateOf(sampler), held)
expect_null(sampler$setState(states$ordered, forceUpdate = TRUE))
expect_identical(cutPointsOf(sampler), attr(held, "cutPoints"))
sampler <- kindOf("whole")
expect_identical(cutPointsOf(sampler)[[1L]], c(0.5, 1.5, 2.5))
expect_identical(withVisible(sampler$setState(states$ordered)), visibly(TRUE))

# a subset rule beside a moved column is as stored: the root splits the
# factor, its left child a threshold the grid lacks with a point left for it
sampler <- dbarts::dbarts(
  y ~ a + b + f,
  frame,
  control = controlWith(n.trees = 5L)
)
state <- stateOf(sampler)
expect_null(attr(state, "cutPoints")[[3L]])
state <- withTree(
  state,
  c(3, 1, -1, -1, -1),
  c(0, fine[30L], 0.1, -0.2, 0.05),
  flags = c(4L, 2L, 0L, 0L, 0L)
)
# levels q and r go right
state[[1L]]$forests[[1L]]$tree.values[1L:8L] <- as.raw(c(6L, rep(0L, 7L)))
sampler$setCutPoints(coarse, 1L)
expect_identical(withVisible(sampler$setState(state)), visibly(FALSE))
expect_null(sampler$setState(state, forceUpdate = TRUE))
moved <- stateOf(sampler)[[1L]]$forests[[1L]]
expect_identical(moved$tree.vars[1L:5L], c(3L, 1L, -1L, -1L, -1L))
expect_identical(moved$tree.flags[1L:5L], as.raw(c(4L, 2L, 0L, 0L, 0L)))
expect_identical(
  moved$tree.values[1L:8L],
  state[[1L]]$forests[[1L]]$tree.values[1L:8L]
)
expect_identical(
  readBin(moved$tree.values[-(1L:8L)], "double", 4L),
  c(coarse[1L], 0.1, -0.2, 0.05)
)

# ---- direction rides the move ---------------------------------------------

# a moved rule keeps the side it sends a missing value to
gone <- seq_len(n) %% 5L == 0L
xHoled <- x
xHoled[gone, 1L] <- NA
holed <- function(...) make(..., data = xHoled)
sampler <- holed(n.trees = 5L, run = FALSE)
onHoled <- cutPointsOf(sampler)[[1L]]
state <- withTree(
  stateOf(sampler),
  c(1, -1, -1),
  c(onHoled[30L], 0.1, -0.2),
  flags = c(3L, 0L, 0L)
)
sampler$setCutPoints(onHoled[c(50L, 75L)], 1L)
expect_identical(withVisible(sampler$setState(state)), visibly(FALSE))
expect_null(sampler$setState(state, forceUpdate = TRUE))
tree <- firstTree(stateOf(sampler))
expect_identical(tree$values, c(onHoled[50L], 0.1, -0.2))
expect_identical(tree$flags, c(3L, 0L, 0L))

# a state stored before the column's first missing value has its sides
# drawn before any split is moved: the draws are the ones a twin makes that
# takes the state on its own grid and has its grid replaced afterwards
donor <- make()
invisible(donor$run(0L, 3L))
stored <- stateOf(donor)
expect_identical(attr(stored, "missing.columns"), c(FALSE, FALSE))
sampler <- holed()
twin <- holed()
sampler$setCutPoints(coarse, 1L)
expect_null(sampler$setState(stored, forceUpdate = TRUE))
twin$setCutPoints(attr(stored, "cutPoints"))
expect_null(twin$setState(stored, forceUpdate = TRUE))
twin$setCutPoints(coarse, 1L, splits = "value")
moved <- stateOf(sampler)
expect_identical(moved, stateOf(twin))
sides <- bitwAnd(as.integer(moved[[1L]]$forests[[1L]]$tree.flags), 1L) == 1L
expect_true(any(sides[moved[[1L]]$forests[[1L]]$tree.vars == 1L]))
expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L))

# ---- linear and gp leaves -------------------------------------------------

# a linear leaf that is not merged keeps its intercept and slope through a
# move
linearWith <- function(...) {
  dbarts::dbarts(
    x,
    y,
    leaf.prior = linear("b"),
    control = controlWith(...)
  )
}
sampler <- linearWith(n.trees = 5L)
state <- withTree(stateOf(sampler), c(1, -1, -1), c(fine[30L], 0.1, -0.2))
slopes <- state[[1L]]$forests[[1L]]$tree.params
expect_identical(length(slopes), 5L)
state[[1L]]$forests[[1L]]$tree.params <- c(0.3, -0.4, slopes[-1L])
sampler$setCutPoints(coarse, 1L)
expect_identical(withVisible(sampler$setState(state)), visibly(FALSE))
expect_null(sampler$setState(state, forceUpdate = TRUE))
moved <- stateOf(sampler)
expect_identical(firstTree(moved)$vars, c(1L, -1L, -1L))
expect_equal(
  firstTree(moved)$values,
  c(coarse[1L], 0.1, -0.2),
  tolerance = 1e-12
)
expect_equal(
  moved[[1L]]$forests[[1L]]$tree.params[1L:2L],
  c(0.3, -0.4),
  tolerance = 1e-12
)

# a gp forest's per-row fits are the stored ones, merged leaves or not
gpWith <- function() {
  sampler <- dbarts::dbarts(
    x,
    y,
    leaf.prior = gp("b"),
    control = controlWith(n.trees = 5L)
  )
  invisible(sampler$run(10L, 1L))
  sampler
}
donor <- gpWith()
stored <- stateOf(donor)
sampler <- gpWith()
invisible(sampler$run(0L, 2L))
sampler$setCutPoints(three, 1L)
expect_identical(withVisible(sampler$setState(stored)), visibly(FALSE))
expect_null(sampler$setState(stored, forceUpdate = TRUE))
moved <- stateOf(sampler)
expect_false(identical(liveOf(moved), liveOf(stored)))
expect_true(all(thresholdsOn(moved, 1L) %in% three))
expect_identical(
  moved[[1L]]$forests[[1L]]$tree.params,
  stored[[1L]]$forests[[1L]]$tree.params
)
expect_true(any(stored[[1L]]$forests[[1L]]$tree.params != 0))
