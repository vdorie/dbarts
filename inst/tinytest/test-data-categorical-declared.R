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
labels.mut <- as.character(g.gap)
labels.mut[1L] <- "d" # the declared-but-unobserved top level
expect_true(sampler.mut$setPredictor(labels.mut, column = 2L))
expect_equal(unname(sampler.mut$data@x[1L, 2L]), 3)
expect_true(all(is.finite(sampler.mut$run(20L, 20L)$train)))

# a column update takes the column's labels, for the test set as for the
# training column, and installs what a whole-frame update installs
sampler.col <- dbarts(train.gap, y.gap, control = control)
sampler.col$setTestPredictor(test.gap)
whole <- sampler.col$data@x.test
sampler.col$setTestPredictor(test.gap[c(2:6, 1L), ])
sampler.col$setTestPredictor(test.gap$g, column = "g")
sampler.col$setTestPredictor(test.gap$x1, column = "x1")
expect_identical(unname(sampler.col$data@x.test), unname(whole))
sampler.col$setTestPredictor(as.character(test.gap$g[c(2:6, 1L)]), column = 2L)
sampler.col$setTestPredictor(as.character(test.gap$g), column = 2L)
expect_identical(unname(sampler.col$data@x.test), unname(whole))
expect_error(
  sampler.col$setTestPredictor(c(0, 1, 2, 3, 3, 0), column = "g"),
  pattern = "column 'g' is categorical; give its values as a factor"
)
expect_error(
  sampler.col$setTestPredictor(c("a", "b", "c", "e", "d", "a"), column = "g"),
  pattern = "column 'g' has label 'e' not among its training levels"
)
expect_error(
  sampler.col$setTestPredictor(c("a", NA, "c", "d", "d", "a"), column = "g"),
  pattern = "column 'g' has missing values, which its training values do not"
)
expect_identical(unname(sampler.col$data@x.test), unname(whole))
sampler.col$setPredictor(
  factor(labels.mut, levels = levels.gap),
  column = "g",
  forceUpdate = TRUE
)
expect_equal(
  unname(sampler.col$data@x[, 2L]),
  match(labels.mut, levels.gap) - 1
)
rm(sampler.col, whole)

# STILL BOUNDED: a label past the declared ones is refused by name, and a
# code at the declared count is out of range on the matrix side
labels.over <- labels.mut
labels.over[1L] <- "e"
expect_error(
  sampler.mut$setPredictor(labels.over, column = 2L),
  pattern = "column 'g' has label 'e' not among its training levels"
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

# a whole-matrix predictor replacement keeps the declared level tables, so a
# sampler whose new rows miss a factor's top level still copies and reloads:
# the level count is the declared one, not the largest code left
local({
  set.seed(3002L)
  n.keep <- 120L
  frame <- data.frame(
    a = runif(n.keep),
    f = factor(sample(letters[1:5], n.keep, TRUE), levels = letters[1:5]),
    o = factor(sample(1:6, n.keep, TRUE), levels = 1:6, ordered = TRUE)
  )
  response <- frame$a +
    as.integer(frame$f) / 3 +
    as.integer(frame$o) / 4 +
    rnorm(n.keep, 0, 0.2)
  sampler <- dbarts(
    response ~ a + f + o,
    frame,
    control = dbartsControl(
      n.trees = 20L,
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 5L,
      n.burn = 50L,
      updateState = FALSE
    )
  )
  invisible(sampler$run())
  levelsBefore <- attr(sampler$data@x, "factor.levels")
  low <- which(as.integer(frame$f) <= 3L & as.integer(frame$o) <= 3L)
  # a matrix of codes is refused on a design with factor columns (dec-B359);
  # the replacement is a data frame, coded by label, and the declared levels
  # stay though the new rows miss the top ones
  replacement <- frame[sample(low, n.keep, TRUE), c("a", "f", "o")]
  expect_error(
    sampler$setPredictor(as.matrix(sampler$data@x)[sample(low, n.keep, TRUE), ]),
    "the predictors 'f', 'o' are factors"
  )
  sampler$setPredictor(replacement, forceUpdate = TRUE)
  expect_identical(attr(sampler$data@x, "factor.levels"), levelsBefore)
  sampler$storeState()
  copied <- sampler$copy()
  expect_true(all(is.finite(copied$run(0L, 2L)$train)))
  reloaded <- unserialize(serialize(sampler, NULL))
  expect_true(all(is.finite(reloaded$run(0L, 2L)$train)))
})
