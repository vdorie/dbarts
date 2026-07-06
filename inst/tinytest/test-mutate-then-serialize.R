# store-after-mutation round trips: every existing state round trip stores a
# never-mutated sampler. Here a warm sampler runs, has a designated column
# replaced (setPredictor forceUpdate), runs again, and stores its state; a cold
# sampler built over the ORIGINAL data then takes the same mutation and the
# stored state and must continue bitwise identically. Covered for linear and gp
# leaves, plus the default constant leaf. n.chains = 2L keeps the draws off R's
# global generator so the continuation is deterministic.

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
  node.prior = linear("x2"),
  control = control
)
invisible(warm.lin$run(20L, 2L))
warm.lin$setPredictor(x2.new, "x2", forceUpdate = TRUE)
invisible(warm.lin$run(0L, 2L))
warm.lin$storeState()

cold.lin <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  node.prior = linear("x2"),
  control = control
)
cold.lin$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.lin$setState(warm.lin$state)
expect_identical(warm.lin$run(0L, 3L), cold.lin$run(0L, 3L))

# gp leaf: same flow
warm.gp <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  node.prior = gp("x2", max.leaf.size = 100L),
  control = control
)
invisible(warm.gp$run(20L, 2L))
warm.gp$setPredictor(x2.new, "x2", forceUpdate = TRUE)
invisible(warm.gp$run(0L, 2L))
warm.gp$storeState()

cold.gp <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  node.prior = gp("x2", max.leaf.size = 100L),
  control = control
)
cold.gp$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.gp$setState(warm.gp$state)
expect_identical(warm.gp$run(0L, 3L), cold.gp$run(0L, 3L))

# default constant leaf: the mutation-then-restore flow with the base prior
warm.const <- dbarts(y ~ x1 + x2 + x3, df, control = control)
invisible(warm.const$run(20L, 2L))
warm.const$setPredictor(x2.new, "x2", forceUpdate = TRUE)
invisible(warm.const$run(0L, 2L))
warm.const$storeState()

cold.const <- dbarts(y ~ x1 + x2 + x3, df, control = control)
cold.const$setPredictor(x2.new, "x2", forceUpdate = TRUE)
cold.const$setState(warm.const$state)
expect_identical(warm.const$run(0L, 3L), cold.const$run(0L, 3L))

rm(
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
