# store-after-mutation round trips: every existing state round trip stores a
# never-mutated sampler. Here a warm sampler runs, has a designated column
# replaced (setPredictor forceUpdate), runs again, and stores its state; a cold
# sampler built over the ORIGINAL data then takes the same mutation and the
# stored state and must reproduce the same model. Covered for linear and gp
# leaves, plus the default constant leaf, and for the copy() and readRDS routes
# that re-create over the mutated data.
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(99)
n <- 150L
x1 <- runif(n)
x2 <- runif(n, -1, 1)
x3 <- runif(n)
mu <- ifelse(x1 > 0.5, x2, 0)
y <- mu + rnorm(n, 0, 0.2)
df <- data.frame(x1, x2, x3, y)
x2.new <- df$x2 * 1.1 # the mutated designated column

control <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L,
  updateState = FALSE
)

# linear leaf: warm-mutate-store, then cold-restore over the original data
warm.lin <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = linear("x2"),
  control = control
)
invisible(warm.lin$run(20L, 2L))
warm.lin$setPredictor(x2.new, "x2", forceUpdate = TRUE)
invisible(warm.lin$run(0L, 2L))
warm.lin$storeState()

cold.lin <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = linear("x2"),
  control = control
)
cold.lin$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.lin$setState(warm.lin$state)
cold.lin$storeState()
statesAgree(cold.lin$state, warm.lin$state)

# gp leaf: same flow
warm.gp <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = gp("x2", max.leaf.size = 50L),
  control = control
)
# the small leaf cap keeps the run fast; the cap's constant fallback is
# reported by each run and pinned here (about half the evaluations at this
# cap, well clear of the quarter that reports it)
expect_warning(
  warm.gp$run(20L, 2L),
  "fell back to a constant leaf",
  strict = TRUE
)
warm.gp$setPredictor(x2.new, "x2", forceUpdate = TRUE)
expect_warning(
  warm.gp$run(0L, 2L),
  "fell back to a constant leaf",
  strict = TRUE
)
warm.gp$storeState()

cold.gp <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = gp("x2", max.leaf.size = 50L),
  control = control
)
cold.gp$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.gp$setState(warm.gp$state)
cold.gp$storeState()
statesAgree(cold.gp$state, warm.gp$state)

# default constant leaf: the mutation-then-restore flow with the base prior
warm.const <- dbarts(y ~ x1 + x2 + x3, df, control = control)
invisible(warm.const$run(20L, 2L))
warm.const$setPredictor(x2.new, "x2", forceUpdate = TRUE)
invisible(warm.const$run(0L, 2L))
warm.const$storeState()

cold.const <- dbarts(y ~ x1 + x2 + x3, df, control = control)
cold.const$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.const$setState(warm.const$state)
cold.const$storeState()
statesAgree(cold.const$state, warm.const$state)

# re-creation over the MUTATED data: copy() and saveRDS/readRDS rebuild from
# data@x, which already holds the replaced column, and must still read the
# saved slopes and kernels through the creation-time standardization the live
# sampler kept, which the state carries
x.new <- data.frame(x1 = runif(20L), x2 = runif(20L, -1, 1), x3 = runif(20L))
recreated <- function(sampler) {
  sampler$storeState()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(sampler, path)
  list(copy = sampler$copy(), reloaded = readRDS(path))
}
keep.control <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 3L,
  keepTrees = TRUE,
  updateState = FALSE
)
for (leaf in c("linear", "gp")) {
  live <- dbarts(
    y ~ x1 + x2 + x3,
    df,
    leaf.prior = if (leaf == "linear") {
      linear("x2")
    } else {
      gp("x2", max.leaf.size = 200L)
    },
    control = keep.control
  )
  invisible(live$run(20L, 3L))
  expect_true(live$setPredictor(2 * x2.new, "x2"), info = leaf)
  invisible(live$run(0L, 3L))
  live.pred <- live$predict(x.new)
  for (route in recreated(live)) {
    expect_identical(route$predict(x.new), live.pred, info = leaf)
  }
  # given the carried rng the continuation is the live one, to the ulps a
  # restore's re-summed fits may differ by
  copied <- live$copy()
  expect_equal(copied$run(0L, 2L)$train, live$run(0L, 2L)$train, info = leaf)
}

# a calibration block the leaf cannot read is refused
live$storeState()
bad <- live$state
bad[[1L]]$forests[[1L]]$leaf.covariate.scale <- 0
expect_error(live$setState(bad), "not consistent with this sampler")
bad <- live$state
bad[[1L]]$forests[[1L]]$leaf.lengthscales <- c(1, 1)
expect_error(live$setState(bad), "not consistent with this sampler")
bad <- live$state
bad[[1L]]$forests[[1L]]$leaf.covariate.center <- "a"
expect_error(live$setState(bad), "'leaf.covariate.center' is malformed")
# a state from before the blocks existed restores from the sampler's own data
old <- live$state
for (chain in seq_along(old)) {
  for (block in c(
    "leaf.covariate.center",
    "leaf.covariate.scale",
    "leaf.lengthscales"
  )) {
    old[[chain]]$forests[[1L]][[block]] <- NULL
  }
}
expect_silent(status <- live$setState(old))
expect_true(status)

rm(
  x.new,
  recreated,
  keep.control,
  leaf,
  live,
  live.pred,
  route,
  copied,
  bad,
  old,
  chain,
  block,
  warm.lin,
  cold.lin,
  warm.gp,
  cold.gp,
  warm.const,
  cold.const,
  control,
  df,
  x1,
  x2,
  x3,
  mu,
  y,
  x2.new,
  n
)
