# A state stored while a column held missing values sends them to a side of
# its rules. Once a forced setPredictor fills them the column routes none, and
# setState, copy(), a reload and a sampler built over the filled rows all
# install the state with those directions dropped, as the forced update
# dropped them from the live trees.
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(83)
n <- 200L
filled <- data.frame(
  x1 = rnorm(n),
  f = factor(sample(letters[1:4], n, replace = TRUE)),
  x2 = rnorm(n)
)
goneX <- seq_len(n) %in% sample(n, 60L)
goneF <- seq_len(n) %in% sample(n, 60L)
y <- ifelse(goneX, 3, filled$x1) +
  ifelse(goneF, -2, as.integer(filled$f) / 2) +
  rnorm(n, sd = 0.2)
holed <- filled
holed$x1[goneX] <- NA
holed$f[goneF] <- NA

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 3L,
  keepTrees = TRUE,
  updateState = FALSE,
  seed = 17L
)
make <- function(data) dbarts(y ~ x1 + f + x2, data, control = control)
fill <- function(sampler) {
  sampler$setPredictor(filled$x1, column = 1L, forceUpdate = TRUE)
  sampler$setPredictor(filled$f, column = 2L, forceUpdate = TRUE)
}
# per live tree of the one forest: its split variables, its split and leaf
# values, and whether each rule sends a missing value right
liveTrees <- function(state) {
  forest <- state[[1L]]$forests[[1L]]
  tree <- rep(seq_along(forest$tree.sizes), forest$tree.sizes)
  list(
    vars = split(forest$tree.vars, tree),
    values = split(forest$tree.values, rep(tree, each = 8L)),
    right = split(bitwAnd(as.integer(forest$tree.flags), 1L) == 1L, tree)
  )
}
treeValues <- function(state) {
  lapply(state, function(chain) chain$forests[[1L]]$tree.values)
}
sendsRight <- function(state) {
  flags <- lapply(state, function(chain) chain$forests[[1L]]$tree.flags)
  bitwAnd(as.integer(unlist(flags)), 1L) == 1L
}
savedFields <- c("saved.vars", "saved.values", "saved.sizes", "saved.flags")

sampler <- make(holed)
invisible(sampler$run(50L, 3L))
sampler$storeState()
stale <- sampler$state
stored <- liveTrees(stale)
# both the ordinal and the categorical column send missing values right
for (column in 1:2) {
  expect_true(any(unlist(stored$right) & unlist(stored$vars) == column))
}
fill(sampler)

# copy() and a reload re-create over the filled predictors from the stale state
duplicate <- sampler$copy()
path <- tempfile(fileext = ".rds")
saveRDS(sampler, path)
reloaded <- readRDS(path)
unlink(path)
other <- make(filled)
invisible(other$run(5L, 1L))

sampler$storeState()
forced <- sampler$state
expect_silent(status <- sampler$setState(stale))
expect_false(status)
expect_silent(status <- other$setState(stale))
expect_false(status)
for (route in list(sampler, duplicate, reloaded, other)) {
  route$storeState()
  statesAgree(route$state, forced)
  expect_identical(treeValues(route$state), treeValues(forced))
  expect_false(any(unlist(liveTrees(route$state)$right)))
  # the saved draws replay by value and keep the directions they were drawn with
  expect_identical(
    route$state[[1L]]$forests[[1L]][savedFields],
    stale[[1L]]$forests[[1L]][savedFields]
  )
  expect_true(all(is.finite(route$run(0L, 3L)$train)))
}

# dropping a direction moves no row: a tree the fill left no empty leaf in is
# the stored tree, rule for rule and leaf for leaf, less its directions
restored <- liveTrees(forced)
whole <- lengths(stored$vars) == lengths(restored$vars)
expect_true(any(vapply(stored$right[whole], any, NA)))
expect_identical(restored$vars[whole], stored$vars[whole])
expect_identical(restored$values[whole], stored$values[whole])

# a column that still holds missing values keeps its directions: the state
# reinstalls as it was stored, and draws the same from it each time
kept <- make(holed)
invisible(kept$run(50L, 3L))
kept$storeState()
held <- kept$state
expect_true(any(unlist(liveTrees(held)$right)))
invisible(kept$run(0L, 3L))
expect_true(kept$setState(held))
kept$storeState()
statesAgree(kept$state, held)
expect_identical(treeValues(kept$state), treeValues(held))
first <- kept$run(0L, 3L)
kept$setState(held)
expect_identical(kept$run(0L, 3L), first)

# two chains, and a sparse-backed column: the stale state restores through
# setState and copy() to the forced update's trees, and both then run
checkRestores <- function(sampler, fill, info) {
  invisible(sampler$run(50L, 3L))
  sampler$storeState()
  stale <- sampler$state
  expect_true(any(sendsRight(stale)), info = info)
  fill(sampler)
  duplicate <- sampler$copy()
  sampler$storeState()
  forced <- sampler$state
  expect_silent(status <- sampler$setState(stale), info = info)
  expect_false(status, info = info)
  for (route in list(sampler, duplicate)) {
    route$storeState()
    statesAgree(route$state, forced)
    expect_identical(treeValues(route$state), treeValues(forced), info = info)
    expect_false(any(sendsRight(route$state)), info = info)
    expect_true(all(is.finite(route$run(0L, 3L)$train)), info = info)
  }
}
twoChains <- control
twoChains@n.chains <- 2L
twoChains@n.threads <- 2L
checkRestores(
  dbarts(y ~ x1 + f + x2, holed, control = twoChains),
  fill,
  "two chains"
)
if (requireNamespace("Matrix", quietly = TRUE)) {
  # mostly zero, so the column stays sparse-backed; its entries carry the holes
  thin <- cbind(ifelse(seq_len(n) %% 4L == 0L, filled$x1, 0), filled$x2)
  thinHoled <- thin
  thinHoled[goneX & thin[, 1L] != 0, 1L] <- NA
  yThin <- ifelse(is.na(thinHoled[, 1L]), 3, thin[, 1L]) + rnorm(n, sd = 0.2)
  checkRestores(
    dbarts(
      Matrix::Matrix(thinHoled, sparse = TRUE),
      yThin,
      control = control,
      sigest = 1
    ),
    function(sampler) {
      sampler$setPredictor(
        Matrix::Matrix(thin, sparse = TRUE),
        forceUpdate = TRUE
      )
    },
    "sparse-backed column"
  )
}
