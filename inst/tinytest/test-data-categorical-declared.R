# declared factor levels are honored end to end: a categorical column's
# category count is the level table its host declares, not the max code any
# training row happens to carry. A factor whose top level goes unobserved in
# training ("a gap factor") therefore keeps a bin for it, and a test or
# mutation row carrying that level is accepted at every entrance - creation,
# setTestPredictor, predict, and setPredictor - exactly as the sparse
# (sparseFactor) route has always accepted it.

set.seed(3001L)
n <- 150L
levels.gap <- c("a", "b", "c", "d")
codes.gap <- sample.int(3L, n, replace = TRUE) # level "d" is never observed
g.gap <- factor(levels.gap[codes.gap], levels = levels.gap)
x1 <- rnorm(n)
y.gap <- 0.6 * codes.gap + x1 + rnorm(n, 0, 0.5)
train.gap <- data.frame(x1 = x1, g = g.gap)

expect_equal(nlevels(g.gap), 4L)
expect_equal(max(as.integer(g.gap)), 3L) # the gap: declared 4, observed 3

control <- dbartsControl(
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
test.gap <- data.frame(
  x1 = rnorm(6L),
  g = factor(c("a", "b", "c", "d", "d", "a"), levels = levels.gap)
)

# CREATION: the design ingests and fits, and the unobserved level owns the
# top code, so getTrees decodes one direction per DECLARED level
sampler.gap <- dbarts(train.gap, y.gap, test = test.gap, control = control)
samples.gap <- sampler.gap$run(20L, 20L)
expect_true(all(is.finite(samples.gap$train)))
expect_true(all(is.finite(samples.gap$test)))

control.keep <- dbartsControl(
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  keepTrees = TRUE,
  updateState = FALSE
)
sampler.keep <- dbarts(train.gap, y.gap, control = control.keep)
invisible(sampler.keep$run(20L, 5L))
trees.gap <- sampler.keep$getTrees()
isRule.g <- trees.gap$var == 2L
expect_true(any(isRule.g))
expect_true(all(nchar(trees.gap$directions[isRule.g]) == 4L))

# SET TEST DATA: a later test set carrying the gap level installs and fits
sampler.set <- dbarts(train.gap, y.gap, control = control)
sampler.set$setTestPredictor(test.gap)
expect_true(all(is.finite(sampler.set$run(20L, 20L)$test)))

# PREDICT: the saved-tree replay routes the gap level too
predictions.gap <- sampler.keep$predict(test.gap)
expect_equal(dim(predictions.gap), c(6L, 5L))
expect_true(all(is.finite(predictions.gap)))

# SET PREDICTOR: installing the gap level into the training column is a valid
# mutation, and the sampler still fits
sampler.mut <- dbarts(train.gap, y.gap, control = control)
codes.mut <- as.double(as.integer(g.gap) - 1L)
codes.mut[1L] <- 3 # the declared-but-unobserved top level
expect_true(sampler.mut$setPredictor(codes.mut, column = 2L))
expect_true(all(is.finite(sampler.mut$run(20L, 20L)$train)))

# STILL BOUNDED: a code at the declared count is out of range on both sides
codes.over <- codes.mut
codes.over[1L] <- 4
expect_error(
  sampler.mut$setPredictor(codes.over, column = 2L),
  pattern = "categorical predictor values must be existing category codes"
)
test.over <- cbind(rnorm(3L), c(0, 1, 4))
colnames(test.over) <- c("x1", "g")
expect_error(
  dbarts(
    dbartsData(train.gap, y.gap, test = test.over),
    control = control
  ),
  pattern = "categorical test predictors must hold existing category codes"
)

# SYMMETRY with the sparse route: the same values given as a sparseFactor
# declare the same K, so both accept the same test set
if (requireNamespace("Matrix", quietly = TRUE)) {
  train.sparse <- data.frame(x1 = x1)
  train.sparse$g <- sparseFactor(g.gap, reference = "b")
  sampler.sparse <- dbarts(
    train.sparse,
    y.gap,
    test = test.gap,
    control = control
  )
  expect_true(all(is.finite(sampler.sparse$run(20L, 20L)$test)))

  # MIXED CONTAINER: a DENSE factor riding beside a sparse column takes its
  # declared count too (the assembleMixedMatrix flavor, which a dense-only
  # design would not exercise)
  train.mixed <- data.frame(x1 = x1, g = g.gap)
  train.mixed$s <- Matrix::sparseVector(
    x = 0.5 + runif(15L),
    i = sort(sample.int(n, 15L)),
    length = n
  )
  data.mixed <- dbartsData(train.mixed, y.gap)
  expect_inherits(data.mixed@x, "dbartsMixedMatrix")
  expect_equal(data.mixed@varTypes, c(0L, 1L, 0L))
  sampler.mixed <- dbarts(data.mixed, control = control)
  expect_true(all(is.finite(sampler.mixed$run(20L, 20L)$train)))

  test.mixed <- data.frame(x1 = rnorm(6L), g = test.gap$g)
  test.mixed$s <- Matrix::sparseVector(x = 1.0, i = 2L, length = 6L)
  sampler.mixed$setTestPredictor(test.mixed)
  expect_true(all(is.finite(sampler.mixed$run(20L, 20L)$test)))
}

# TIER CROSSING: the declared count decides the 63/64 inline/pooled boundary,
# so a design whose declared table crosses it while its observed codes do not
# pools its rule masks. A state saved before declared levels were honored -
# reproduced here by a sampler over the same values under the trimmed level
# table, which stays inline - is refused rather than silently restored.
set.seed(3002L)
n.wide <- 300L
levels.wide <- sprintf("V%02d", 1:70)
codes.wide <- sample.int(60L, n.wide, replace = TRUE) # 61..70 unobserved
g.declared <- factor(levels.wide[codes.wide], levels = levels.wide)
g.trimmed <- factor(levels.wide[codes.wide], levels = levels.wide[1:60])
z.wide <- rnorm(n.wide)
y.wide <- ifelse(codes.wide >= 30L, 1.5, 0) + z.wide + rnorm(n.wide, 0, 0.4)

control.wide <- dbartsControl(
  n.trees = 15L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  keepTrees = TRUE,
  updateState = TRUE
)
sampler.declared <- dbarts(
  data.frame(z = z.wide, g = g.declared),
  y.wide,
  control = control.wide
)
invisible(sampler.declared$run(20L, 5L))
sampler.trimmed <- dbarts(
  data.frame(z = z.wide, g = g.trimmed),
  y.wide,
  control = control.wide
)
invisible(sampler.trimmed$run(20L, 5L))

# the declared design pools (a mask side channel), the trimmed one stays inline
expect_true(length(sampler.declared$state[[1L]]$forests[[1L]]$tree.masks) > 0L)
expect_null(sampler.trimmed$state[[1L]]$forests[[1L]]$tree.masks)

expect_error(
  sampler.declared$setState(sampler.trimmed$state),
  pattern = "missing required block 'tree.masks'"
)

# SAME TIER: a state saved by another declared-table sampler restores
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
state.declared <- sampler.declared$state
sampler.restore <- dbarts(
  data.frame(z = z.wide, g = g.declared),
  y.wide,
  control = control.wide
)
sampler.restore$setState(state.declared)
sampler.restore$storeState()
statesAgree(sampler.restore$state, state.declared)

rm(
  sampler.gap,
  sampler.keep,
  sampler.set,
  sampler.mut,
  sampler.declared,
  sampler.trimmed,
  sampler.restore,
  samples.gap,
  trees.gap,
  isRule.g,
  predictions.gap,
  state.declared,
  control,
  control.keep,
  control.wide,
  train.gap,
  test.gap,
  test.over,
  codes.gap,
  codes.mut,
  codes.over,
  codes.wide,
  levels.gap,
  levels.wide,
  g.gap,
  g.declared,
  g.trimmed,
  x1,
  y.gap,
  y.wide,
  z.wide,
  n,
  n.wide
)
