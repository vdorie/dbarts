# A state stored before a forced setPredictor can route no row of the new
# predictors to some of its leaves. setState, copy(), a reload and a same-grid
# warm start all install it, merging those leaves into their parents exactly as
# the forced update merged them, as 0.9-34's restore did.
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(71)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L)
y <- 4 * (x[, 1L] > 0.5) + x[, 2L] + rnorm(n, sd = 0.2)
z <- rbinom(n, 1L, 0.5)
labels <- factor(sample(0:2, n, replace = TRUE))
# the extremes stay, so the cut grid does not move; the interior collapses
xNew <- x
interior <- x[, 1L] > min(x[, 1L]) & x[, 1L] < max(x[, 1L])
xNew[interior, 1L] <- 0.5

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 3L,
  keepTrees = TRUE,
  updateState = FALSE,
  seed = 17L
)
# the merged leaf values, which statesAgree leaves out
treeValues <- function(state) {
  lapply(state, function(chain) lapply(chain$forests, `[[`, "tree.values"))
}
noEmptyLeaf <- function(sampler) {
  trees <- sampler$getTrees()
  !any(trees$var == -1L & trees$n == 0L)
}

checkStaleRestore <- function(make, info) {
  sampler <- make(x)
  invisible(sampler$run(40L, 3L))
  sampler$storeState()
  stale <- sampler$state
  sampler$setPredictor(xNew, forceUpdate = TRUE)

  # copy() and a reload re-create over the new predictors from the stale state
  duplicate <- sampler$copy()
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  reloaded <- readRDS(path)
  unlink(path)

  sampler$storeState()
  forced <- sampler$state
  expect_false(
    statesAgree(forced, stale, expect = FALSE),
    info = paste(info, "the forced update merged a leaf")
  )
  expect_silent(sampler$setState(stale))
  sampler$storeState()
  statesAgree(sampler$state, forced)
  expect_identical(treeValues(sampler$state), treeValues(forced), info = info)
  for (route in list(duplicate, reloaded)) {
    route$storeState()
    statesAgree(route$state, forced)
    expect_identical(treeValues(route$state), treeValues(forced), info = info)
  }
  for (route in list(sampler, duplicate, reloaded)) {
    result <- route$run(0L, 3L)
    expect_true(all(is.finite(result$train)), info = info)
    expect_true(noEmptyLeaf(route), info = info)
  }
}

checkStaleRestore(
  function(x) dbarts(x, y, control = control),
  "gaussian:"
)
checkStaleRestore(
  function(x) dbarts(x, labels, family = "multinomial", control = control),
  "multinomial:"
)
checkStaleRestore(
  function(x) {
    dbarts(
      x,
      y,
      forests = list(forest(), forest(basis = ~ factor(z))),
      control = control
    )
  },
  "bcf:"
)

# a restore with every leaf occupied merges nothing: it reinstalls the stored
# trees as they were
exact <- dbarts(x, y, control = control)
invisible(exact$run(40L, 3L))
exact$storeState()
held <- exact$state
invisible(exact$run(0L, 3L))
exact$setState(held)
exact$storeState()
statesAgree(exact$state, held)

# a merge weighs a subtree's leaves under the state's own latents, not the
# destination's: a sampler whose latents moved after the forced update and a
# copy() starting cold restore the same stale state to the same trees, leaf
# values and next draws. Whether a merged subtree's leaf weights differ is
# draw-dependent, so each family runs over several seeds
yCount <- rpois(n, exp(x[, 1L] + x[, 2L]))
makeLatentSampler <- function(name, seed) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 3L,
    updateState = FALSE,
    seed = seed
  )
  switch(
    name,
    logistic = dbarts(
      x,
      as.integer(y > 2.5),
      family = binomial(link = "logit"),
      control = control
    ),
    nbinom = dbarts(x, yCount, family = "nbinom", control = control),
    student = dbarts(x, y, family = student(df = 5), control = control)
  )
}
for (name in c("logistic", "nbinom", "student")) {
  for (seed in 1:6) {
    info <- paste(name, seed)
    moved <- makeLatentSampler(name, seed)
    invisible(moved$run(40L, 3L))
    moved$storeState()
    stale <- moved$state
    moved$setPredictor(xNew, forceUpdate = TRUE)
    invisible(moved$run(0L, 5L))
    cold <- moved$copy()
    moved$setState(stale)
    moved$storeState()
    cold$storeState()
    statesAgree(moved$state, cold$state)
    expect_identical(
      treeValues(moved$state),
      treeValues(cold$state),
      info = info
    )
    expect_identical(moved$run(0L, 3L), cold$run(0L, 3L), info = info)
  }
}

# a same-grid warm start whose donor leaves a leaf with no rows: the
# destination's trees are merged, so it carries no empty leaf, and it restores
# and copies itself
donor <- dbarts(x, y, control = control)
invisible(donor$run(40L, 3L))
donor$storeState()
destination <- dbarts(xNew, y, control = control)
destination$installTrees(donor)
invisible(destination$run(0L, 1L))
expect_true(noEmptyLeaf(destination))
destination$storeState()
expect_silent(destination$setState(destination$state))
destinationCopy <- destination$copy()
expect_true(all(is.finite(destinationCopy$run(0L, 3L)$train)))
