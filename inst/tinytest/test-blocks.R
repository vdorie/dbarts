# Block-additive constraints (variant A, docs/design/interaction-constraints.md):
# the blocks() surface and its fit-time validation, a constrained fit confining
# every tree to one declared group, and the per-forest BCF split. The engine
# confinement + warm-start refusal gates live in tests/cpp/test_state.cpp
# (testBlockAdditiveConfinement).

set.seed(22L)

# group a getTrees data.frame into one tree per (chain, sample, tree)
splitTrees <- function(trees) {
  cols <- intersect(c("chain", "sample", "tree"), names(trees))
  split(trees, trees[cols], drop = TRUE)
}
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
doFit <- function(blocks = NULL) {
  do.call(bart2, c(list(y ~ x1 + x2 + x3, df, blocks = blocks), fitArgs))
}

# ---- a blocks() fit runs and confines every tree to its group -----------------

# groups {x1} and {x2, x3} -> columns {1} and {2, 3}
groups <- list(1L, c(2L, 3L))
fit <- doFit(blocks(groups = list("x1", c("x2", "x3"))))
trees <- extract(fit, type = "trees")
expect_true(someSplit(trees))
expect_true(allConfined(trees, groups))

# an explicit trees.per.group runs and stays confined
fitSplit <- doFit(blocks(
  groups = list("x1", c("x2", "x3")),
  trees.per.group = c(5L, 15L)
))
expect_true(allConfined(extract(fitSplit, type = "trees"), groups))

# a numeric-index partition is equivalent
fitIdx <- doFit(blocks(groups = list(1L, c(2L, 3L))))
expect_true(allConfined(extract(fitIdx, type = "trees"), groups))

# ---- fit-time validation (safe over fast) -------------------------------------

expect_error(blocks(), "requires 'groups'")

# an un-named predictor would be masked out of every tree and go dead
expect_error(
  doFit(blocks(groups = list("x1", "x2"))),
  "name every predictor"
)
# a predictor named in two groups
expect_error(
  doFit(blocks(groups = list(c("x1", "x2"), c("x2", "x3")))),
  "disjoint"
)
# an unrecognized name
expect_error(
  doFit(blocks(groups = list("x1", c("x2", "x3"), "nope"))),
  "unrecognized variable name 'nope'"
)
# trees.per.group of the wrong length
expect_error(
  doFit(blocks(groups = list("x1", c("x2", "x3")), trees.per.group = 20L)),
  "one entry per group"
)
# trees.per.group that does not sum to n.trees
expect_error(
  doFit(blocks(
    groups = list("x1", c("x2", "x3")),
    trees.per.group = c(5L, 5L)
  )),
  "sum to"
)
# a non-positive capacity
expect_error(
  doFit(blocks(
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
invisible(dbarts:::bartcoreRun(bc, 150L, 0L))
tauTrees <- dbarts:::bartcoreGetTrees(
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
