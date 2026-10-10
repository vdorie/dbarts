# setState's two forms. Without forceUpdate a state that fits the sampler is
# installed and TRUE is returned, visibly; a state one of whose trees would
# have to be changed to fit is not installed, FALSE is returned and the
# sampler is left as it was. With forceUpdate = TRUE the state is installed,
# repaired where it must be, and NULL is returned invisibly. copy() and a
# reload force, silently. The kept draws are output: a state installs whatever
# store it was taken from, the newest of its draws that fit.

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(4127L, sample.kind = "Rejection")
n <- 150L
columns <- c("a", "b", "c")
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, columns))
y <- 3 * (x[, 1L] > 0.5) + x[, 2L] + sin(4 * x[, 3L]) + rnorm(n, sd = 0.2)
z <- rep_len(c(0, 1), n)
xTest <- matrix(runif(30L), 10L, 3L, dimnames = list(NULL, columns))

controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 4L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = 31L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
# a sampler past its burn-in holding its own stored state; built twice the
# same way, two are the same sampler
warmed <- function(..., control = controlWith(), data = x, response = y) {
  sampler <- dbarts::dbarts(data, response, ..., control = control)
  if (control@keepTrees) {
    invisible(sampler$run(20L, control@n.samples))
  } else {
    invisible(sampler$run(20L, 0L))
  }
  sampler$storeState()
  sampler
}
visibly <- function(value) list(value = value, visible = TRUE)
forcedValue <- list(value = NULL, visible = FALSE)
keptDraws <- function(sampler) sampler$predict(xTest)
liveTrees <- function(sampler) sampler$getTrees(current = TRUE)
isLive <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_isValidPointer, sampler$pointer)
}

# A state that goes in without force: TRUE, visibly, and NULL invisibly with
# it, the two samplers then drawing the same 20 draws. Returns the first.
expectClean <- function(make, state, info) {
  unforced <- make()
  forced <- make()
  expect_identical(withVisible(unforced$setState(state)), visibly(TRUE), info)
  expect_identical(unforced$state, state, info = info)
  expect_identical(
    withVisible(forced$setState(state, forceUpdate = TRUE)),
    forcedValue,
    info = info
  )
  expect_identical(liveTrees(unforced), liveTrees(forced), info = info)
  expect_identical(unforced$run(0L, 20L), forced$run(0L, 20L), info = info)
  unforced
}
# A state declined without force: FALSE, visibly, and the sampler is the twin
# that made no call - its state field, its trees and kept draws, the whole
# state its engine then stores (sigma, k, cut points, latents and generator
# among it) and its next 20 draws. Forced, the state goes in, NULL invisibly.
# Returns the forced sampler, its state stored.
expectDeclined <- function(make, state, info) {
  sampler <- make()
  twin <- make()
  forced <- make()
  field <- sampler$state
  expect_identical(withVisible(sampler$setState(state)), visibly(FALSE), info)
  expect_identical(sampler$state, field, info = info)
  expect_identical(liveTrees(sampler), liveTrees(twin), info = info)
  if (sampler$control@keepTrees) {
    expect_identical(keptDraws(sampler), keptDraws(twin), info = info)
  }
  sampler$storeState()
  twin$storeState()
  expect_identical(sampler$state, twin$state, info = info)
  expect_identical(sampler$run(0L, 20L), twin$run(0L, 20L), info = info)
  expect_identical(
    withVisible(forced$setState(state, forceUpdate = TRUE)),
    forcedValue,
    info = info
  )
  expect_identical(forced$state, state, info = info)
  forced$storeState()
  forced
}

# ---- the argument ---------------------------------------------------------

own <- warmed()
for (bad in list(NA, "partial", c(TRUE, FALSE), 1, NULL)) {
  expect_error(
    own$setState(own$state, forceUpdate = bad),
    pattern = "'forceUpdate' must be TRUE or FALSE",
    fixed = TRUE
  )
}

# ---- clean ----------------------------------------------------------------

# the sampler's own state, and a recorded chain put into a fresh sampler
expectClean(warmed, own$state, "own state")
expectClean(
  function() dbarts::dbarts(x, y, control = controlWith(seed = 77L)),
  own$state,
  "a recorded chain in a fresh sampler"
)

# a state stored before a predictor change that empties no leaf
nudged <- x[, 2L] + runif(n, -1e-3, 1e-3)
expectClean(
  function() {
    sampler <- warmed()
    expect_true(sampler$setPredictor(nudged, 2L))
    sampler
  },
  own$state,
  "after a predictor change that empties no leaf"
)

# latents redrawn under other weights (Student-t) and other censoring (aft)
weights <- rep_len(c(1, 4, 4), n)
studentWith <- function(w) {
  function() {
    sampler <- warmed(
      weights = weights,
      family = dbarts:::student(5),
      control = controlWith(keepTrees = FALSE)
    )
    sampler$setWeights(w)
    sampler
  }
}
expectClean(
  studentWith(replace(weights, 1:20, 0)),
  studentWith(weights)()$state,
  "other weights"
)
time <- exp(y / 3)
bound <- as.double(quantile(time, 0.7))
aftWith <- function(status) {
  function() {
    warmed(
      response = cbind(pmin(time, bound), status),
      family = "aft",
      control = controlWith(keepTrees = FALSE)
    )
  }
}
expectClean(
  aftWith(rep(1, n)),
  aftWith(as.double(time <= bound))()$state,
  "other censoring"
)

# a column filled after its state was stored keeps its rules' sides for
# missing values: the sampler has seen it missing, and no row moves leaf
xMissing <- x
xMissing[seq(3L, n, by = 9L), 2L] <- NA
withMissing <- warmed(data = xMissing)
expect_true(any(withMissing$state[[1L]]$forests[[1L]]$tree.flags == as.raw(3L)))
expectClean(
  function() {
    sampler <- warmed(data = xMissing)
    expect_true(sampler$setPredictor(x[, 2L], 2L))
    sampler
  },
  withMissing$state,
  "a filled column's directions"
)

# ---- the kept store -------------------------------------------------------

# a ring of ten that has wrapped: thirteen draws leave the oldest at slot 4
ten <- warmed(control = controlWith(n.samples = 10L))
invisible(ten$run(0L, 3L))
ten$storeState()
expect_identical(attr(ten$state, "currentSampleNum"), 3L)
tenDraws <- keptDraws(ten)

# ten into four: the newest four, in order; n.samples does not move
expectClean(warmed, ten$state, "ten kept draws into four")
four <- warmed()
expect_true(four$setState(ten$state))
expect_identical(four$control@n.samples, 4L)
expect_identical(keptDraws(four), tenDraws[, 7:10])
invisible(four$run(0L, 1L))
expect_identical(keptDraws(four)[, 1:3], tenDraws[, 8:10])

# four, wrapped, into ten: all four in order, and three sweeps make seven
small <- warmed()
invisible(small$run(0L, 2L))
small$storeState()
fourDraws <- keptDraws(small)
expectClean(
  function() warmed(control = controlWith(n.samples = 10L)),
  small$state,
  "four kept draws into ten"
)
large <- warmed(control = controlWith(n.samples = 10L))
expect_true(large$setState(small$state))
expect_identical(keptDraws(large), fourDraws)
invisible(large$run(0L, 3L))
expect_identical(dim(keptDraws(large)), c(10L, 7L))
expect_identical(keptDraws(large)[, 1:4], fourDraws)

# four into none, and none into four
bare <- function() warmed(control = controlWith(keepTrees = FALSE))
expectClean(bare, small$state, "four kept draws into none")
expectClean(warmed, bare()$state, "no kept draws into four")
emptied <- warmed()
expect_true(emptied$setState(bare()$state))
expect_error(keptDraws(emptied), pattern = "holds no recorded draws")
invisible(emptied$run(0L, 2L))
expect_identical(dim(keptDraws(emptied)), c(10L, 2L))

# an equal store continues the ring where the state left it: a copy and its
# source, five sweeps on, hold the same draws in the same slots and store the
# same state, to the rounding a re-creation's rebuilt fits leave
duplicate <- ten$copy()
expect_equal(duplicate$run(0L, 5L), ten$run(0L, 5L), tolerance = 1e-10)
expect_equal(keptDraws(duplicate), keptDraws(ten), tolerance = 1e-10)
duplicate$storeState()
ten$storeState()
statesAgree(duplicate$state, ten$state)
savedShape <- function(state) {
  c(
    state[[1L]]$forests[[1L]][c("saved.vars", "saved.sizes", "saved.flags")],
    attributes(state)[c("currentSampleNum", "recordedDraws")]
  )
}
expect_identical(savedShape(duplicate$state), savedShape(ten$state))
expect_identical(attr(ten$state, "currentSampleNum"), 8L)

# ---- not clean: a leaf no row reaches -------------------------------------

# the extremes stay, so the cut grid does not move; the interior collapses
xCollapsed <- x
interior <- x[, 1L] > min(x[, 1L]) & x[, 1L] < max(x[, 1L])
xCollapsed[interior, 1L] <- 0.5
afterCollapse <- function(...) {
  function() {
    sampler <- warmed(...)
    sampler$setPredictor(xCollapsed, forceUpdate = TRUE)
    sampler
  }
}
# forced, the install merges as the forced update merged
merged <- afterCollapse()()
merged$storeState()
forced <- expectDeclined(afterCollapse(), own$state, "a leaf no row reaches")
expect_false(identical(
  merged$state[[1L]]$forests[[1L]]$tree.sizes,
  own$state[[1L]]$forests[[1L]]$tree.sizes
))
expect_identical(
  forced$state[[1L]]$forests[[1L]][c("tree.vars", "tree.sizes", "tree.values")],
  merged$state[[1L]]$forests[[1L]][c("tree.vars", "tree.sizes", "tree.values")]
)

# in the variance forest alone: the forced update's state with the stored
# variance trees put back
varianceFields <- paste0("variance.", c("vars", "values", "sizes", "flags"))
heteroscedastic <- afterCollapse(variance = TRUE)
merged <- heteroscedastic()
merged$storeState()
staleVariance <- merged$state
staleVariance[[1L]][varianceFields] <- warmed(variance = TRUE)$state[[
  1L
]][varianceFields]
expect_false(identical(
  staleVariance[[1L]]$variance.sizes,
  merged$state[[1L]]$variance.sizes
))
forced <- expectDeclined(
  heteroscedastic,
  staleVariance,
  "a variance leaf no row reaches"
)
expect_identical(
  forced$state[[1L]][varianceFields],
  merged$state[[1L]][varianceFields]
)

# a dead pointer stays dead after a state is declined, and the field stays
path <- tempfile(fileext = ".rds")
saveRDS(afterCollapse()(), path)
reloaded <- readRDS(path)
unlink(path)
expect_false(isLive(reloaded))
expect_identical(reloaded$state, own$state)
expect_identical(reloaded$setState(own$state), FALSE)
expect_false(isLive(reloaded))
expect_identical(reloaded$state, own$state)
expect_null(reloaded$setState(own$state, forceUpdate = TRUE))
expect_true(isLive(reloaded))

# ---- not clean: a side for missing values on a column never missing -------

# another sampler's state, whose column held missing values: this sampler's
# never has, so the sides its rules record cannot stand
forced <- expectDeclined(
  warmed,
  withMissing$state,
  "a direction on a column never missing"
)
expect_false(any(forced$state[[1L]]$forests[[1L]]$tree.flags == as.raw(3L)))

# ---- not clean, and repaired: a tree its forest does not allow ------------

# Tree 1 of a stored block, split variables and values (cut points and leaf
# values as stored), and the block with a hand-built tree in its place.
firstTree <- function(block, prefix) {
  size <- block[[paste0(prefix, ".sizes")]][1L]
  list(
    vars = block[[paste0(prefix, ".vars")]][seq_len(size)],
    values = readBin(block[[paste0(prefix, ".values")]], "double", size)
  )
}
withFirstTree <- function(block, prefix, tree) {
  field <- function(name) paste0(prefix, ".", name)
  old <- seq_len(block[[field("sizes")]][1L])
  block[[field("vars")]] <- c(tree$vars, block[[field("vars")]][-old])
  block[[field("values")]] <- c(
    writeBin(tree$values, raw()),
    block[[field("values")]][-seq_len(length(old) * 8L)]
  )
  # 2 tags a threshold split
  block[[field("flags")]] <- c(
    as.raw(2L * (tree$vars > 0L)),
    block[[field("flags")]][-old]
  )
  block[[field("sizes")]][1L] <- length(tree$vars)
  block
}
cutOf <- function(j) {
  cuts <- attr(own$state, "cutPoints")[[j]]
  cuts[length(cuts) %/% 2L]
}
cut <- vapply(1:3, cutOf, 0)
left <- x[, 1L] <= cut[1L]
weightedMean <- function(values, rows) sum(values * rows) / sum(rows)
# the mean tree: a at the root; on its left c over a leaf and a split on b;
# a leaf on its right
leaves <- c(0.03, -0.05, 0.11, -0.02)
reach <- c(
  sum(left & x[, 3L] <= cut[3L]),
  sum(left & x[, 3L] > cut[3L] & x[, 2L] <= cut[2L]),
  sum(left & x[, 3L] > cut[3L] & x[, 2L] > cut[2L])
)
expect_true(all(reach > 0L) && length(unique(reach)) == 3L)
meanTree <- list(
  vars = c(1L, 3L, -1L, 2L, -1L, -1L, -1L),
  values = c(cut[1L], cut[3L], leaves[1L], cut[2L], leaves[2:4])
)
belowRoot <- list(
  vars = c(1L, -1L, -1L),
  values = c(cut[1L], weightedMean(leaves[1:3], reach), leaves[4L])
)
belowSecond <- list(
  vars = c(1L, 3L, -1L, -1L, -1L),
  values = c(
    cut[1L],
    cut[3L],
    leaves[1L],
    weightedMean(leaves[2:3], reach[2:3]),
    leaves[4L]
  )
)
# the variance tree: a at the root, b on its left, the merged leaf at the
# geometric mean
scales <- c(0.7, 1.9, 1.1)
reachVariance <- c(
  sum(left & x[, 2L] <= cut[2L]),
  sum(left & x[, 2L] > cut[2L])
)
varianceTree <- list(
  vars = c(1L, 2L, -1L, -1L, -1L),
  values = c(cut[1L], cut[2L], scales)
)
varianceBelowRoot <- list(
  vars = c(1L, -1L, -1L),
  values = c(
    cut[1L],
    exp(weightedMean(log(scales[1:2]), reachVariance)),
    scales[3L]
  )
)

forest <- dbarts::dbartsForests$forest
blocks <- dbarts::dbartsForests$blocks
interactions <- dbarts::dbartsForests$interactions
onAB <- c("a", "b")
cases <- list(
  list(
    info = "a split off the forest's columns",
    make = function() warmed(forests = list(forest(vars = onAB))),
    forest = 1L,
    expected = belowRoot
  ),
  list(
    info = "a split off a moderator forest's columns",
    # its kept draws have no combined prediction to compare
    make = function() {
      warmed(
        forests = list(forest(), forest(basis = z, vars = onAB)),
        control = controlWith(keepTrees = FALSE)
      )
    },
    forest = 2L,
    expected = belowRoot
  ),
  list(
    info = "a split off a tree's block",
    make = function() {
      groups <- blocks(groups = list("a", c("b", "c")))
      warmed(forests = list(forest(blocks = groups)))
    },
    forest = 1L,
    expected = belowRoot
  ),
  list(
    info = "a second variable under max.order = 1",
    make = function() warmed(interactions = interactions(max.order = 1L)),
    forest = 1L,
    expected = belowRoot
  ),
  list(
    info = "a forbidden pair on one path",
    make = function() {
      warmed(interactions = interactions(forbid = list(c(1L, 2L))))
    },
    forest = 1L,
    expected = belowSecond
  ),
  list(
    info = "a variance split off the variance forest's columns",
    make = function() warmed(variance = "a"),
    forest = 0L,
    expected = varianceBelowRoot
  )
)
# the hand-built tree in a sampler's stored state, and the tree a state holds
# in its place
handState <- function(state, forest) {
  if (forest == 0L) {
    state[[1L]] <- withFirstTree(state[[1L]], "variance", varianceTree)
  } else {
    state[[1L]]$forests[[forest]] <- withFirstTree(
      state[[1L]]$forests[[forest]],
      "tree",
      meanTree
    )
  }
  state
}
treeIn <- function(state, forest) {
  if (forest == 0L) {
    firstTree(state[[1L]], "variance")
  } else {
    firstTree(state[[1L]]$forests[[forest]], "tree")
  }
}
expectTree <- function(got, expected, info) {
  expect_identical(got$vars, expected$vars, info = info)
  expect_equal(got$values, expected$values, tolerance = 1e-12, info = info)
}
for (case in cases) {
  state <- handState(case$make()$state, case$forest)
  # collapsed at the first split from the root the forest does not allow:
  # the split above it kept, everything beneath it one leaf at the
  # row-weighted mean, and every other tree as stored
  forced <- expectDeclined(case$make, state, case$info)
  expectTree(treeIn(forced$state, case$forest), case$expected, case$info)
  restored <- handState(forced$state, case$forest)
  expect_identical(restored[[1L]], state[[1L]], info = case$info)
  expect_true(forced$setState(forced$state), info = case$info)
  expect_true(all(is.finite(forced$run(0L, 5L)$train)), info = case$info)

  # a copy and a reload of a sampler holding such a state repair it, silently
  holder <- case$make()
  holder$state <- state
  expect_silent(duplicate <- holder$copy())
  path <- tempfile(fileext = ".rds")
  saveRDS(holder, path)
  reloaded <- readRDS(path)
  unlink(path)
  for (route in list(duplicate, reloaded)) {
    expect_silent(route$storeState())
    expectTree(treeIn(route$state, case$forest), case$expected, case$info)
  }
}
# where nothing forbids them the hand-built trees are clean
unrestricted <- function() warmed(variance = TRUE)
state <- handState(handState(unrestricted()$state, 1L), 0L)
expectClean(unrestricted, state, "the hand-built trees, unrestricted")
taken <- unrestricted()
expect_true(taken$setState(state))
taken$storeState()
expectTree(treeIn(taken$state, 1L), meanTree, "unrestricted")
expectTree(treeIn(taken$state, 0L), varianceTree, "unrestricted")

# a monotone tree whose leaf values, as stored, fall where they must rise:
# reseeded, every leaf 0
monotone <- function() warmed(monotone = c(a = "increasing"))
falling <- list(vars = c(1L, -1L, -1L), values = c(cut[1L], 0.2, -0.2))
state <- monotone()$state
state[[1L]]$forests[[1L]] <- withFirstTree(
  state[[1L]]$forests[[1L]],
  "tree",
  falling
)
forced <- expectDeclined(monotone, state, "a monotone tree out of order")
expectTree(
  treeIn(forced$state, 1L),
  list(vars = falling$vars, values = c(cut[1L], 0, 0)),
  "a monotone tree out of order"
)

# ---- the kept store of a factor of more than 63 levels --------------------

# Such a factor's rules keep their levels in words beside each tree, kept
# draws included: a store of another size takes the words with the draws.
wide <- data.frame(a = x[, 1L], f = factor(rep_len(sprintf("l%02d", 1:70), n)))
yWide <- 3 * (wide$a > 0.5) + 2 * (as.integer(wide$f) <= 30L) + sin(seq_len(n))
wideTest <- wide[seq(1L, n, by = 7L), ]
warmedWide <- function(n.samples, more = 0L) {
  sampler <- dbarts::dbarts(
    yWide ~ a + f,
    wide,
    control = controlWith(n.samples = n.samples)
  )
  invisible(sampler$run(40L, n.samples))
  if (more > 0L) {
    invisible(sampler$run(0L, more))
  }
  sampler$storeState()
  sampler
}
# the bytes of the words the kept draws in `slots` of a stored state hold,
# in that order: ten trees to a slot, two words of eight bytes to a rule on f
keptWords <- function(state, slots) {
  forest <- state[[1L]]$forests[[1L]]
  slot <- rep((seq_along(forest$saved.sizes) - 1L) %/% 10L, forest$saved.sizes)
  owner <- rep(slot[forest$saved.vars == 2L], each = 16L)
  unlist(lapply(slots, function(s) forest$saved.masks[owner == s]))
}
# the kept draws of `from` put into `into`: the words of the draws that
# moved, slots `newest` of the source, and then what the draws predict, a
# rule read against other words reading past them
expectWordsMoved <- function(into, from, newest, draws, info) {
  expect_true(into$setState(from$state), info = info)
  into$storeState()
  moved <- keptWords(into$state, seq_along(newest) - 1L)
  expected <- keptWords(from$state, newest)
  expect_true(length(moved) > 0L, info = info)
  expect_identical(moved, expected, info = info)
  if (!identical(moved, expected)) {
    return(invisible(NULL))
  }
  expect_identical(
    into$predict(wideTest),
    from$predict(wideTest)[, draws],
    info = info
  )
  trees <- from$getTrees()
  trees <- trees[trees$sample %in% draws, ]
  trees$sample <- trees$sample - min(draws) + 1L
  expect_identical(as.list(into$getTrees()), as.list(trees), info = info)
}
# ten that have wrapped into four: the oldest of the newest four is at slot 9
tenWide <- warmedWide(10L, 3L)
expect_identical(attr(tenWide$state, "currentSampleNum"), 3L)
expectWordsMoved(warmedWide(4L), tenWide, c(9L, 0:2), 7:10, "ten into four")
# four that have wrapped into ten: the oldest is at slot 2
fourWide <- warmedWide(4L, 2L)
expect_identical(attr(fourWide$state, "currentSampleNum"), 2L)
expectWordsMoved(warmedWide(10L), fourWide, c(2:3, 0:1), 1:4, "four into ten")

# ---- a decline reconciles nothing -----------------------------------------

# A state stored under other weights or another censoring has its latent
# values redrawn once it is installed. A declined state is not installed: a
# redraw would move the latent values and the generators of a sampler the
# call says it left alone.
binary <- as.double(y > median(y))
logisticUnder <- function(w, collapse) {
  function() {
    sampler <- warmed(
      response = binary,
      weights = w,
      family = binomial(link = "logit"),
      control = controlWith(keepTrees = FALSE)
    )
    if (collapse) {
      sampler$setWeights(rev(w))
      sampler$setPredictor(xCollapsed, forceUpdate = TRUE)
    }
    sampler
  }
}
aftCollapsed <- function(status) {
  function() {
    sampler <- aftWith(status)()
    sampler$setPredictor(xCollapsed, forceUpdate = TRUE)
    sampler
  }
}
cases <- list(
  list(
    info = "declined under other weights",
    make = logisticUnder(weights, TRUE),
    state = logisticUnder(weights, FALSE)()$state,
    digest = "weights.digest"
  ),
  list(
    info = "declined under other censoring",
    # the sampler that holds censored rows: the redraw is of their times
    make = aftCollapsed(as.double(time <= bound)),
    state = aftWith(rep(1, n))()$state,
    digest = "survival.digest"
  )
)
for (case in cases) {
  held <- case$make()
  held$storeState()
  expect_false(
    identical(attr(held$state, case$digest), attr(case$state, case$digest)),
    info = case$info
  )
  expectDeclined(case$make, case$state, case$info)
}

# ---- declined on another cut grid -----------------------------------------

# a state is judged on its own cut points, which are in place while it is;
# declined, the sampler's own come back with everything else
coarse <- afterCollapse(control = controlWith(n.cuts = 20L))
expect_false(identical(
  attr(coarse()$state, "cutPoints"),
  attr(own$state, "cutPoints")
))
forced <- expectDeclined(coarse, own$state, "a state on another cut grid")
expect_identical(attr(forced$state, "cutPoints"), attr(own$state, "cutPoints"))
