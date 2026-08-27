# Block-additive constraints (variant A): the blocks() surface and its
# fit-time validation, a constrained fit confining
# every tree to one declared group, and the per-forest BCF split. The engine
# confinement + warm-start refusal gates live in tests/cpp/test_state.cpp
# (testBlockAdditiveConfinement).

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(22L)

# group a getTrees data.frame into one tree per (chain, sample, tree)
source(
  system.file("common", "splitTrees.R", package = "dbarts"),
  local = TRUE
)
# the distinct 1-based split variables a flat tree uses (empty for a stump)
usedVars <- function(k) {
  v <- k$var
  sort(unique(v[!is.na(v) & v > 0]))
}
# every tree's variable set must be a subset of exactly one declared group; a
# never-split tree (empty set) is a subset of all and passes
allConfined <- function(trees, groups) {
  all(vapply(
    splitTrees(trees),
    function(k) {
      used <- usedVars(k)
      length(used) == 0L ||
        any(vapply(groups, function(g) all(used %in% g), logical(1)))
    },
    logical(1)
  ))
}
# does at least one tree actually split (so confinement is non-trivial)?
someSplit <- function(trees) {
  any(vapply(
    splitTrees(trees),
    function(k) length(usedVars(k)) > 0L,
    logical(1)
  ))
}

# both blocks carry signal, so trees in each group have something to split on
n <- 200L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- 2 * x1 + ifelse(x2 > 0.5, 1, -1) + 0.7 * x3 + rnorm(n, 0, 0.3)
df <- data.frame(y, x1, x2, x3)

fitArgs <- list(
  n.trees = 20L,
  n.samples = 15L,
  n.burn = 40L,
  n.chains = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
doFitBlocks <- function(blocks = NULL) {
  do.call(bart2, c(list(y ~ x1 + x2 + x3, df, blocks = blocks), fitArgs))
}

# ---- a blocks() fit runs and confines every tree to its group -----------------

# groups {x1} and {x2, x3} -> columns {1} and {2, 3}
groups <- list(1L, c(2L, 3L))
fit <- doFitBlocks(blocks(groups = list("x1", c("x2", "x3"))))
trees <- extract(fit, type = "trees")
expect_true(someSplit(trees))
expect_true(allConfined(trees, groups))

# an explicit trees.per.group runs and stays confined
fitSplit <- doFitBlocks(blocks(
  groups = list("x1", c("x2", "x3")),
  trees.per.group = c(5L, 15L)
))
expect_true(allConfined(extract(fitSplit, type = "trees"), groups))

# a numeric-index partition is equivalent
fitIdx <- doFitBlocks(blocks(groups = list(1L, c(2L, 3L))))
expect_true(allConfined(extract(fitIdx, type = "trees"), groups))

# ---- a never-split tree is a subset of all blocks, not a violation -------

# x4 carries no signal, so at least one tree confined to its own block stays a
# stump (root leaf, no split) after burn-in; allConfined must accept the empty
# variable set rather than reject it
x4 <- runif(n)
dfC3 <- data.frame(y, x1, x2, x3, x4)
groupsC3 <- list(1L, c(2L, 3L), 4L)
fitC3 <- do.call(
  bart2,
  c(
    list(
      y ~ x1 + x2 + x3 + x4,
      dfC3,
      blocks = blocks(
        groups = list("x1", c("x2", "x3"), "x4"),
        trees.per.group = c(7L, 7L, 7L)
      )
    ),
    modifyList(fitArgs, list(n.trees = 21L))
  )
)
treesC3 <- extract(fitC3, type = "trees")
byTreeC3 <- splitTrees(treesC3)
x4TreeIdx <- vapply(byTreeC3, function(k) k$tree[1], integer(1)) %in% 15:21
expect_true(any(vapply(
  byTreeC3[x4TreeIdx],
  function(k) length(usedVars(k)) == 0L,
  logical(1)
)))
expect_true(allConfined(treesC3, groupsC3))

# ---- per-block tree counts match trees.per.group -------------------------------

# the deterministic contiguous assignment: the first counts[1] trees belong
# to group 1, the next counts[2] to group 2, and so on. A tree that splits must
# stay within ITS OWN assigned block - a stronger check than allConfined's
# "any block" test.
blockOfTreeIndex <- function(idx, counts) {
  findInterval(idx - 1L, cumsum(counts)) + 1L
}
checkFixedCapacity <- function(trees, groups, counts) {
  byTree <- splitTrees(trees)
  idx <- vapply(byTree, function(k) k$tree[1], integer(1))
  blk <- blockOfTreeIndex(idx, counts)
  all(mapply(
    function(k, b) {
      used <- usedVars(k)
      length(used) == 0L || all(used %in% groups[[b]])
    },
    byTree,
    blk
  ))
}

# default even split of 20 trees over 2 groups: 10 and 10
expect_true(checkFixedCapacity(
  extract(fit, type = "trees"),
  groups,
  c(10L, 10L)
))
# the explicit trees.per.group from above: 5 and 15
expect_true(checkFixedCapacity(
  extract(fitSplit, type = "trees"),
  groups,
  c(5L, 15L)
))

# ---- fit-time validation (safe over fast) -------------------------------------

expect_error(blocks(), "requires 'groups'")

# an un-named predictor would be masked out of every tree and go dead
expect_error(
  doFitBlocks(blocks(groups = list("x1", "x2"))),
  "name every predictor"
)
# a predictor named in two groups
expect_error(
  doFitBlocks(blocks(groups = list(c("x1", "x2"), c("x2", "x3")))),
  "disjoint"
)
# an unrecognized name
expect_error(
  doFitBlocks(blocks(groups = list("x1", c("x2", "x3"), "nope"))),
  "unrecognized variable name 'nope'"
)
# trees.per.group of the wrong length
expect_error(
  doFitBlocks(blocks(
    groups = list("x1", c("x2", "x3")),
    trees.per.group = 20L
  )),
  "one entry per group"
)
# trees.per.group that does not sum to n.trees
expect_error(
  doFitBlocks(blocks(
    groups = list("x1", c("x2", "x3")),
    trees.per.group = c(5L, 5L)
  )),
  "sum to"
)
# a non-positive capacity
expect_error(
  doFitBlocks(blocks(
    groups = list("x1", c("x2", "x3")),
    trees.per.group = c(0L, 20L)
  )),
  "positive"
)
# more groups than trees, with the default even split
expect_error(
  do.call(
    bart2,
    c(
      list(
        y ~ x1 + x2 + x3,
        df,
        blocks = blocks(groups = list("x1", "x2", "x3"))
      ),
      modifyList(fitArgs, list(n.trees = 2L))
    )
  ),
  "at least one tree"
)

# ---- blocks() and interactions() compose on one forest ------------------------

fitBoth <- do.call(
  bart2,
  c(
    list(
      y ~ x1 + x2 + x3,
      df,
      blocks = blocks(groups = list("x1", c("x2", "x3"))),
      interactions = interactions(max.order = 1)
    ),
    fitArgs
  )
)
expect_true(allConfined(extract(fitBoth, type = "trees"), groups))

# ---- blocks() and monotone() compose on one forest -----------------------------

# monotone constrains leaf-value direction, blocks() constrains split
# selection - orthogonal seams; the fit must run and confinement must still
# hold
fitMono <- do.call(
  bart2,
  c(
    list(
      y ~ x1 + x2 + x3,
      df,
      blocks = blocks(groups = list("x1", c("x2", "x3"))),
      monotone = c(x1 = "+", x3 = "+")
    ),
    fitArgs
  )
)
expect_true(allConfined(extract(fitMono, type = "trees"), groups))

# ---- a fixed blocks() config + seed reproduces draws ---------------------------

set.seed(77L)
fitRepA <- doFitBlocks(blocks(groups = list("x1", c("x2", "x3"))))
set.seed(77L)
fitRepB <- doFitBlocks(blocks(groups = list("x1", c("x2", "x3"))))
expect_identical(
  extract(fitRepA, type = "trees"),
  extract(fitRepB, type = "trees")
)
expect_identical(fitRepA$sigma, fitRepB$sigma)
expect_identical(fitRepA$yhat.train, fitRepB$yhat.train)

# ---- per-forest BCF: tau confined to a moderator block ------------------------

p <- 4L
xb <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
z <- rbinom(n, 1L, 0.5)
tau <- ifelse(xb[, 1] > 0.5, 2, -1) + ifelse(xb[, 2] > 0.5, 1, -1)
yb <- 2 * sin(pi * xb[, 3]) + xb[, 4] + z * tau + rnorm(n, sd = 0.2)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE
)
bcSampler <- dbarts(xb, yb, control = control)
bc <- dbarts:::bartcoreBCFSampler(
  bcSampler,
  z,
  n.trees.treatment = 20L,
  moderators = c(1L, 2L),
  tau.blocks = blocks(groups = list("x1", "x2"))
)
invisible(bartcoreRun(bc, 150L, 0L))
tauTrees <- bartcoreGetTrees(
  bc,
  chainNums = 1L,
  treeNums = 1:20,
  current = TRUE,
  forest = 1L
)
# every tau tree confined to {x1} or {x2}
expect_true(allConfined(tauTrees, list(1L, 2L)))

# tau.blocks must partition exactly the moderator set: naming a non-moderator errors
targetSampler <- dbarts(xb, yb, control = control)
expect_error(
  dbarts:::bartcoreBCFSampler(
    targetSampler,
    z,
    n.trees.treatment = 20L,
    moderators = c(1L, 2L),
    tau.blocks = blocks(groups = list("x1", c("x2", "x3")))
  ),
  "available predictors"
)

# ---- warm-start refusal + pre-mutation integrity -------------------------------

# an unconstrained donor is free to mix x1 with x2/x3 in any tree; a target
# whose block 0 confines 15 of its 20 trees to {x1} alone almost certainly
# receives a donor tree that splits outside that mask, and
# columnMaskStateFeasible (chain.hpp) must refuse the whole install before
# touching any live state - the same mechanism and error text as the BCF
# columnMask refusal test in test-interactions.R.
wsControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.burn = 40L,
  n.samples = 15L,
  updateState = TRUE,
  keepTrees = FALSE
)
wsDonor <- dbarts(y ~ x1 + x2 + x3, df, control = wsControl)
invisible(wsDonor$run())

wsTarget <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  control = wsControl,
  blocks = blocks(
    groups = list("x1", c("x2", "x3")),
    trees.per.group = c(15L, 5L)
  )
)
expect_error(wsTarget$installTrees(wsDonor), "column restriction")

# the refusal happened before any live state mutation: the target sampler is
# still usable afterwards
wsSamples <- wsTarget$run()
expect_equal(dim(wsSamples$train), c(n, 15L))
