source(
  system.file("common", "leafPriorChecks.R", package = "dbarts"),
  local = TRUE
)

library(dbarts, quietly = TRUE)

set.seed(99)
n <- 250L
x1 <- runif(n)
x2 <- runif(n, -1, 1)
x3 <- runif(n)
g <- factor(sample(letters[1:3], n, replace = TRUE))
mu <- ifelse(x1 > 0.5, x2, 0)
y <- mu + rnorm(n, 0, 0.2)
df <- data.frame(x1, x2, x3, y)

# designation validation happens when the prior resolves against the data
expect_error(
  dbarts(y ~ x1 + x2 + x3, df, leaf.prior = linear("zz")),
  pattern = "unrecognized column"
)
expect_error(
  dbarts(y ~ x1 + x2 + x3, df, leaf.prior = linear(10)),
  pattern = "out of range"
)
expect_error(
  dbarts(y ~ x1 + x2 + x3, df, leaf.prior = linear(c(2, 2))),
  pattern = "duplicate"
)
# a fractional numeric columns index is refused, naming the argument, rather
# than silently truncated (coerceOrError's integer branch)
expect_error(
  dbarts(y ~ x1 + x2 + x3, df, leaf.prior = linear(1.5)),
  "'columns' must be a whole number; got '1.5'",
  fixed = TRUE
)
expect_error(
  dbarts(y ~ x1 + x2 + g, data.frame(df, g), leaf.prior = linear("g")),
  pattern = "must be continuous"
)
# an unresolved designation cannot enter a model object directly
expect_error(
  new("dbartsModel", leaf.prior = dbarts:::linear("x2")),
  pattern = "resolved against data"
)

# fitting recovers the varying slope structure
control <- dbartsControl(
  n.trees = 20L,
  n.chains = 1L,
  n.samples = 40L,
  n.burn = 150L,
  keepTrees = TRUE,
  updateState = FALSE
)
set.seed(0)
sampler <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  test = df[1:5, c("x1", "x2", "x3")],
  leaf.prior = linear("x2"),
  control = control
)
samples <- sampler$run()
fits <- rowMeans(samples$train)
expect_true(sum((fits - mu)^2) < 0.2 * sum((mean(y) - mu)^2))

# recorded test fits match a saved-tree replay of the same rows
predictions <- sampler$predict(as.matrix(df[1:5, c("x1", "x2", "x3")]))
expect_equal(predictions, samples$test, tolerance = 1e-10)

# getTrees reports one slope column per covariate, NA on internal nodes
trees <- sampler$getTrees(treeNums = 1:2, sampleNums = 1L)
expect_true("beta.x2" %in% names(trees))
expect_true(all(is.na(trees$beta.x2[trees$var > 0])))
expect_true(all(!is.na(trees$beta.x2[trees$var == -1])))

# plotTree labels linear leaves with their coefficients; the leaf
# covariate designation is fixed at creation: a replacement model with a
# constant leaf prior is refused
checkPlotTreeAndFixedPrior(sampler)

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
# state serialization carries the slope arrays: a restored sampler
# reproduces the model
list2env(
  # unlike a literal leaf.prior = linear(...) written directly inside a
  # dbarts() call, a linear() object threaded through this helper's own
  # 'leaf.prior' parameter does not reach dbarts()'s NSE routing, so it
  # must be built through the qualified name here
  checkStateRoundTrip(
    y ~ x1 + x2 + x3,
    df,
    dbarts:::linear("x2"),
    10L
  ),
  environment()
)

# the mutable-data surface stays live under linear leaves: a mutated
# sampler's continued fit agrees with a from-scratch fit of the mutated
# data (setPredictor re-quantizes differently from a fresh build, so the
# comparison is statistical, not identical), and the live trees stay
# well-formed
set.seed(1)
sampler.mut <- dbarts(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = linear("x2"),
  control = control
)
invisible(sampler.mut$run(50L, 5L))
x2.new <- df$x2 * 1.1
expect_silent(sampler.mut$setPredictor(x2.new, "x2", forceUpdate = TRUE))
more <- sampler.mut$run(0L, 200L)
fits.mut <- rowMeans(more$train)
liveTrees <- sampler.mut$getTrees(current = TRUE)
expect_true(
  all(is.na(liveTrees$beta.x2[liveTrees$var > 0])) &&
    all(!is.na(liveTrees$beta.x2[liveTrees$var == -1]))
)

df.mut <- df
df.mut$x2 <- x2.new
set.seed(101)
sampler.fresh <- dbarts(
  y ~ x1 + x2 + x3,
  df.mut,
  leaf.prior = linear("x2"),
  control = control
)
fits.fresh <- rowMeans(sampler.fresh$run(150L, 200L)$train)
# rmse between independently-seeded fits of the same mutated data, 60
# seeds: mean 0.041, sd 0.0076; bound clears mean + 4 sd (0.072)
expect_true(sqrt(mean((fits.mut - fits.fresh)^2)) < 0.08)

# a probit response composes with linear leaves
z <- rbinom(n, 1L, pnorm(mu / 0.5))
df.binary <- data.frame(x1, x2, x3, z)
set.seed(2)
sampler.binary <- dbarts(
  z ~ x1 + x2 + x3,
  df.binary,
  leaf.prior = linear("x2"),
  control = control
)
samples.binary <- sampler.binary$run(100L, 20L)
expect_true(all(is.finite(samples.binary$train)))

# linear leaves ride the data-handle views: a full-rows view matches the
# raw-data path bitwise, standardizing with the parent's constants; a
# proper fold serves its held-out rows through the gathered covariates
# threaded through checkDataHandleViews()'s own 'leaf.prior' parameter, so
# linear() must be the qualified name (see the checkStateRoundTrip() call
# above)
list2env(
  checkDataHandleViews(
    y ~ x1 + x2 + x3,
    df,
    dbarts:::linear("x2"),
    15L,
    n,
    mu
  ),
  environment()
)

# views still refuse raw-predictor mutation under linear leaves
x.view <- as.matrix(sampler.view$data@x)
storage.mode(x.view) <- "double"
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_setPredictor,
    fold$ptr,
    x.view,
    FALSE,
    0L
  ),
  pattern = "views hold none"
)

# xbart accepts a linear leaf prior, with its k standing in for a missing
# k argument and the k grid overriding per cell
xbart.linear <- xbart(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = linear("x2", k = 3),
  n.samples = 60L,
  n.burn = c(60L, 30L),
  n.reps = 2L,
  n.trees = 15L,
  n.threads = 1L,
  seed = 1L
)
expect_true(all(is.finite(xbart.linear)))
xbart.grid <- xbart(
  y ~ x1 + x2 + x3,
  df,
  leaf.prior = linear("x2"),
  k = c(1, 4),
  n.samples = 60L,
  n.burn = c(60L, 30L),
  n.reps = 2L,
  n.trees = 15L,
  n.threads = 1L,
  seed = 1L,
  drop = FALSE
)
expect_equal(dim(xbart.grid)[3L], 2L)
expect_error(
  xbart(
    y ~ x1 + x2 + g,
    data.frame(df, g),
    leaf.prior = linear("g"),
    n.threads = 1L
  ),
  pattern = "must be continuous"
)

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
  df,
  df.mut,
  df.binary,
  x1,
  x2,
  x3,
  g,
  y,
  z,
  mu,
  x2.new,
  n,
  control.view,
  sampler.view,
  handle,
  view,
  full,
  samples.view,
  samples.full,
  testRows,
  fold,
  x.view,
  samples.fold,
  xbart.linear,
  xbart.grid
)

# a refused whole-matrix replacement leaves the leaf covariates as they were:
# a constant x1 empties every leaf a split on it bounds, so the change rolls
# back, and the sampler then continues as a twin that never attempted it.
# The rollback repartitions, so the fits agree to rounding, not bitwise.
set.seed(711)
n <- 120L
df.roll <- data.frame(x1 = runif(n), x2 = runif(n, -1, 1))
df.roll$y <- ifelse(df.roll$x1 > 0.5, df.roll$x2, 0) + rnorm(n, 0, 0.2)
control.roll <- dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 20L,
  updateState = FALSE
)
makeRollSampler <- function() {
  set.seed(712)
  sampler <- dbarts(
    y ~ x1 + x2,
    df.roll,
    leaf.prior = linear("x2"),
    control = control.roll
  )
  set.seed(713)
  invisible(sampler$run(30L, 0L))
  sampler
}
refused <- makeRollSampler()
twin <- makeRollSampler()
expect_false(refused$setPredictor(
  cbind(x1 = rep(0.5, n), x2 = -df.roll$x2),
  forceUpdate = FALSE
))
set.seed(714)
samples.refused <- refused$run(0L, 20L)
set.seed(714)
samples.twin <- twin$run(0L, 20L)
expect_equal(samples.refused$train, samples.twin$train, tolerance = 1e-10)
rm(n, df.roll, control.roll, makeRollSampler, refused, twin)
rm(samples.refused, samples.twin)

# no cap on leaf covariates: 9 and 12 columns fit, and the saved-tree replay
# at the training rows is the recorded train fits
set.seed(41)
xWide <- matrix(runif(200L * 12L), 200L, 12L)
colnames(xWide) <- paste0("w", 1:12)
yWide <- xWide[, 1L] + xWide[, 2L] * (xWide[, 3L] > 0.5) + rnorm(200L, 0, 0.2)
for (q in c(9L, 12L)) {
  set.seed(1)
  wideSampler <- dbarts(
    xWide,
    yWide,
    leaf.prior = linear(seq_len(q)),
    control = dbartsControl(
      n.trees = 10L,
      n.chains = 1L,
      n.samples = 20L,
      n.burn = 20L,
      keepTrees = TRUE,
      updateState = FALSE
    )
  )
  wideSamples <- wideSampler$run()
  expect_true(all(is.finite(wideSamples$train)))
  expect_true(all(is.finite(wideSamples$sigma)))
  expect_equal(
    wideSampler$predict(xWide),
    wideSamples$train,
    tolerance = 1e-12
  )
}
