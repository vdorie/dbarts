# A zero split probability excludes its column as a column list does, on
# every path: the forest never splits on it, $setModel refuses a zero on a
# column a live tree splits on, and a state or donor that splits on it is
# declined unforced and cut back forced.

forest <- dbarts::dbartsForests$forest
cgm <- dbarts::dbartsPriors$cgm

set.seed(41L)
n <- 150L
x <- cbind(a = rnorm(n), b = rnorm(n), c = rnorm(n))
y <- x[, "a"] + 2 * x[, "c"] + rnorm(n, 0, 0.5)

makeControl <- function(n.chains = 1L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 1L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 13L
  )
}
# sigest is fixed: the residual scale estimate does not differ by prior here,
# but the comparisons below are stated at one estimate
build <- function(..., design = x, n.chains = 1L) {
  dbarts(design, y, ..., control = makeControl(n.chains), sigest = 1)
}
zeroLast <- cgm(split.probs = c(0.5, 0.5, 0))
splits <- function(draws) rowSums(draws$varcount)
splitsOnLast <- function(sampler, current = FALSE) {
  sum(sampler$getTrees(current = current)$var == 3L)
}

# the zero draws what the column list draws, and the stated probabilities are
# reported as written
zeroed <- build(tree.prior = zeroLast)
listed <- build(tree.prior = zeroLast, forests = list(forest(vars = 1:2)))
zeroDraws <- zeroed$run(0L, 200L)
listDraws <- listed$run(0L, 200L)
expect_identical(zeroDraws$train, listDraws$train)
expect_identical(zeroDraws$sigma, listDraws$sigma)
expect_identical(zeroDraws$varcount, listDraws$varcount)
expect_identical(unname(splits(zeroDraws)[3L]), 0)
expect_true(all(splits(zeroDraws)[1:2] > 0))
expect_identical(zeroed$model@tree.prior@splitProbabilities, c(0.5, 0.5, 0))

# before any run $setModel takes a zero and draws what creating with it draws
installed <- build()
installed$setModel(build(tree.prior = zeroLast)$model)
expect_identical(
  installed$run(0L, 50L)$train,
  build(tree.prior = zeroLast)$run(0L, 50L)$train
)

# a zero on a column some chain's live trees split on is refused, naming it,
# and the sampler is left as a twin never asked; by position where the
# columns have no names
refusedZero <- function(design, pattern) {
  sampler <- build(design = design, n.chains = 2L)
  twin <- build(design = design, n.chains = 2L)
  sampler$run(0L, 30L)
  twin$run(0L, 30L)
  expect_true(splitsOnLast(sampler, current = TRUE) > 0)
  before <- sampler$model
  expect_error(
    sampler$setModel(build(design = design, tree.prior = zeroLast)$model),
    pattern
  )
  expect_identical(sampler$model, before)
  expect_identical(sampler$run(0L, 10L)$train, twin$run(0L, 10L)$train)
}
refusedZero(x, "split probability zero to the predictor \"c\": a live tree")
refusedZero(unname(x), "split probability zero to the predictor in column 3")

# a state or donor that splits on the column is declined, and the sampler is
# left as it was; forced, every such split is cut back
live <- build(tree.prior = zeroLast)
invisible(live$run(0L, 50L))
donor <- build()
invisible(donor$run(0L, 50L))
live$storeState()
donor$storeState()
expect_true(splitsOnLast(donor, current = TRUE) > 0)
stateBefore <- live$state
treesBefore <- live$getTrees()
expect_identical(live$setState(donor$state), FALSE)
expect_identical(live$installTrees(donor), FALSE)
expect_identical(live$state, stateBefore)
expect_identical(live$getTrees(), treesBefore)
forced <- live$copy()
expect_null(forced$setState(donor$state, forceUpdate = TRUE))
expect_identical(splitsOnLast(forced, current = TRUE), 0L)
expect_true(any(forced$getTrees(current = TRUE)$var > 0L))
expect_identical(unname(splits(forced$run(0L, 50L))[3L]), 0)
warm <- live$copy()
expect_null(warm$installTrees(donor, forceUpdate = TRUE))
expect_identical(splitsOnLast(warm, current = TRUE), 0L)
expect_true(any(warm$getTrees(current = TRUE)$var > 0L))
expect_identical(unname(splits(warm$run(0L, 50L))[3L]), 0)
