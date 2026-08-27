source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "leafPriorChecks.R", package = "dbarts"),
  local = TRUE
)

library(dbarts, quietly = TRUE)

set.seed(99)
n <- 250L
x1 <- runif(n)
x2 <- runif(n, -1, 1)
g <- factor(sample(letters[1:3], n, replace = TRUE))
mu <- sin(4 * pi * x1)
y <- mu + rnorm(n, 0, 0.2)
df <- data.frame(x1, x2, y)
z <- rbinom(n, 1L, pnorm(2 * mu / 3))
df.binary <- data.frame(x1, x2, z)

# designation and argument validation happens when the prior resolves
expect_error(
  dbarts(y ~ x1 + x2, df, node.prior = gp("zz", k = 2)),
  pattern = "unrecognized column"
)
expect_error(
  dbarts:::gp("x1", k = 2, lengthscale = -1),
  pattern = "must be positive"
)
expect_error(
  dbarts:::gp("x1", k = 2, max.leaf.size = 0L),
  pattern = "positive integer"
)
expect_error(
  dbarts(y ~ x1 + x2, df, node.prior = gp("x1", k = 2, lengthscale = c(1, 2))),
  pattern = "length 1 or match"
)
expect_error(
  dbarts(y ~ x1 + x2 + g, data.frame(df, g), node.prior = gp("g", k = 2)),
  pattern = "must be continuous"
)
expect_error(
  new("dbartsModel", node.prior = dbarts:::gp("x1")),
  pattern = "resolved against data"
)

# the chi hyperprior composes with gp leaves: k is sampled from the drawn
# functions' standardized magnitudes
control.chi <- dbartsControl(n.trees = 10L, n.chains = 1L, updateState = FALSE)
set.seed(3)
sampler.chi <- dbarts(
  y ~ x1 + x2,
  df,
  node.prior = gp("x1", k = chi(1.25), max.leaf.size = 100L),
  control = control.chi
)
samples.chi <- sampler.chi$run(100L, 20L)
expect_true(all(is.finite(samples.chi$k)) && all(samples.chi$k > 0))
expect_true(length(unique(samples.chi$k)) > 1L)

# binary responses default k to chi under gp leaves like everywhere else
set.seed(4)
sampler.chi.binary <- dbarts(
  z ~ x1 + x2,
  df.binary,
  node.prior = gp("x1", max.leaf.size = 100L),
  control = control.chi
)
samples.chi.binary <- sampler.chi.binary$run(60L, 10L)
expect_true(all(is.finite(samples.chi.binary$train)))
expect_true(length(unique(samples.chi.binary$k)) > 1L)

# zero weights are ignored like the other leaf types: no likelihood
# contribution, fits ride the conditional and stay finite
w0 <- rep_len(c(1, 2, 1, 1, 0), n)
set.seed(6)
sampler.w0 <- suppressWarnings(
  dbarts(
    y ~ x1 + x2,
    df,
    weights = w0,
    node.prior = gp("x1", k = 2, max.leaf.size = 100L),
    control = control.chi
  )
)
samples.w0 <- sampler.w0$run(60L, 10L)
expect_true(all(is.finite(samples.w0$train)))
expect_true(all(is.finite(samples.w0$sigma)))

# fitting recovers the smooth signal; the small cap keeps large leaves on
# the constant fallback while splits grow
control <- dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.samples = 40L,
  n.burn = 150L,
  keepTrees = TRUE,
  updateState = FALSE
)
set.seed(0)
sampler <- dbarts(
  y ~ x1 + x2,
  df,
  test = df[1:5, c("x1", "x2")],
  node.prior = gp("x1", max.leaf.size = 100L),
  control = control
)
samples <- sampler$run()
fits <- rowMeans(samples$train)
expect_true(sum((fits - mu)^2) < 0.2 * sum((mean(y) - mu)^2))

# recorded test fits match a saved-tree replay of the same rows: the saved
# blocks hold the exact weights the live evaluation used
predictions <- sampler$predict(as.matrix(df[1:5, c("x1", "x2")]))
expect_equal(predictions, samples$test, tolerance = 1e-12)

# getTrees reports no coefficient columns and NA leaf values: the function
# rides prediction only
trees <- sampler$getTrees(treeNums = 1:2, sampleNums = 1L)
expect_false(any(startsWith(names(trees), "beta.")))
expect_true(all(is.na(trees$value[trees$var == -1L])))
expect_true(all(!is.na(trees$value[trees$var > 0L])))

# plotTree labels gp leaves; the leaf covariate designation and kind are
# fixed at creation
checkPlotTreeAndFixedPrior(sampler)

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
# state serialization carries the fits slabs and saved blocks: a restored
# sampler reproduces the model
list2env(
  # unlike a literal node.prior = gp(...) written directly inside a dbarts()
  # call, a gp() object threaded through this helper's own 'node.prior'
  # parameter does not reach dbarts()'s NSE routing, so it must be built
  # through the qualified name here
  checkStateRoundTrip(
    y ~ x1 + x2,
    df,
    dbarts:::gp("x1", max.leaf.size = 100L),
    5L
  ),
  environment()
)

# the mutable-data surface stays live under gp leaves, including the
# designated column itself: a mutated sampler's continued fit agrees with
# a from-scratch fit of the mutated data (setPredictor re-quantizes
# differently from a fresh build, so the comparison is statistical, not
# identical), and the live trees stay well-formed
set.seed(1)
sampler.mut <- dbarts(
  y ~ x1 + x2,
  df,
  node.prior = gp("x1", max.leaf.size = 100L),
  control = control
)
invisible(sampler.mut$run(50L, 5L))
x1.new <- pmin(df$x1 * 1.1, 1)
expect_silent(sampler.mut$setPredictor(x1.new, "x1", forceUpdate = TRUE))
more <- sampler.mut$run(0L, 200L)
fits.mut <- rowMeans(more$train)
liveTrees <- sampler.mut$getTrees(current = TRUE)
expect_true(
  all(is.na(liveTrees$value[liveTrees$var == -1L])) &&
    all(!is.na(liveTrees$value[liveTrees$var > 0L]))
)

df.mut <- df
df.mut$x1 <- x1.new
set.seed(101)
sampler.fresh <- dbarts(
  y ~ x1 + x2,
  df.mut,
  node.prior = gp("x1", max.leaf.size = 100L),
  control = control
)
fits.fresh <- rowMeans(sampler.fresh$run(150L, 200L)$train)
# rmse between independently-seeded fits of the same mutated data, 30
# seeds: mean 0.084, sd 0.0093; bound clears mean + 4 sd (0.121)
expect_true(sqrt(mean((fits.mut - fits.fresh)^2)) < 0.15)

# a probit response composes with gp leaves under an explicit fixed k
set.seed(2)
sampler.binary <- dbarts(
  z ~ x1 + x2,
  df.binary,
  node.prior = gp("x1", k = 2, max.leaf.size = 100L),
  control = control
)
samples.binary <- sampler.binary$run(100L, 20L)
expect_true(all(is.finite(samples.binary$train)))

# gp leaves ride the data-handle views: a full-rows view matches the
# raw-data path bitwise, standardizing with the parent's constants; a
# proper fold serves its held-out rows through the gathered covariates
# threaded through checkDataHandleViews()'s own 'node.prior' parameter, so
# gp() must be the qualified name (see the checkStateRoundTrip() call above)
list2env(
  checkDataHandleViews(
    y ~ x1 + x2,
    df,
    dbarts:::gp("x1", max.leaf.size = 100L),
    10L,
    n,
    mu
  ),
  environment()
)

# xbart accepts a gp node prior, with its k standing in for a missing k
# argument
xbart.gp <- xbart(
  y ~ x1 + x2,
  df,
  node.prior = gp("x1", k = 3, max.leaf.size = 100L),
  n.samples = 60L,
  n.burn = c(60L, 30L, 0L),
  n.reps = 2L,
  n.trees = 10L,
  n.threads = 1L,
  seed = 1L
)
expect_true(all(is.finite(xbart.gp)))

# sparse predictor matrices hold no raw covariate values
if (requireNamespace("Matrix", quietly = TRUE)) {
  x.sparse <- Matrix::Matrix(as.matrix(df[, c("x1", "x2")]), sparse = TRUE)
  expect_error(
    dbarts(x.sparse, y, node.prior = gp(1L, k = 2)),
    pattern = "sparse predictor matrices"
  )
}

rm(
  sampler,
  sampler.state,
  sampler.restored,
  sampler.mut,
  sampler.fresh,
  sampler.binary,
  samples,
  samples.binary,
  more,
  fits.mut,
  fits.fresh,
  liveTrees,
  trees,
  predictions,
  fits,
  control,
  control.state,
  control.view,
  sampler.view,
  handle,
  view,
  full,
  samples.view,
  samples.full,
  testRows,
  fold,
  samples.fold,
  xbart.gp,
  df,
  df.mut,
  df.binary,
  x1,
  x2,
  g,
  y,
  z,
  mu,
  x1.new,
  n,
  control.chi,
  sampler.chi,
  samples.chi,
  sampler.chi.binary,
  samples.chi.binary,
  w0,
  sampler.w0,
  samples.w0
)
