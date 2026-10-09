# A copy installs the stored state as a reload does: a state whose saved-tree
# store has another capacity than the control now names (n.samples edited
# after the state was stored) is taken with its own capacity, not refused.

set.seed(0)
n <- 80L
x <- matrix(rnorm(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n)

make <- function() {
  dbarts::dbarts(
    x,
    y,
    control = dbarts::dbartsControl(
      keepTrees = TRUE,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 10L,
      n.samples = 20L,
      n.burn = 5L,
      verbose = FALSE,
      seed = 3L
    )
  )
}

# the control slot edited in place, the state stored before it
sampler <- make()
sampler$run(2L, 5L)
sampler$control@n.samples <- 8L
dupe <- sampler$copy()
expect_true(inherits(dupe, "dbartsSampler"))
expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)

# the supported route, setControl, after the store
sampler <- make()
sampler$run(2L, 5L)
newControl <- sampler$control
newControl@n.samples <- 8L
sampler$setControl(newControl)
dupe <- sampler$copy()
expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)

# an untouched seeded sampler without a stored state starts as its source does
fresh <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    updateState = FALSE,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 20L,
    verbose = FALSE,
    seed = 3L
  )
)
dupe <- fresh$copy()
expect_equal(dupe$run(3L, 5L)$train, fresh$run(3L, 5L)$train)
