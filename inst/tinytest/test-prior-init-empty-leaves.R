# sampleTreesFromPrior draws the tree prior CONDITIONED on the empty-leaf-free
# set the move kernels price, by per-tree
# rejection. The conditioning predicate is the veto's own - membership, not
# weight - so every leaf of a from-prior forest holds a row, a leaf may hold
# only zero-weight rows, and the draw does not read the weights at all: a
# sampler carrying a zero-weight half-space and one carrying no weights draw
# the same forests from the same seed.
#
# The oracle needs no tree walk: the drawn trees report, per node, the rows
# that reach it (getTrees) and, routing ONLY the positive-weight rows through
# them (getTrees(newdata = )), the rows the likelihood can see.

set.seed(20260818L)
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
w <- as.double(x[, 1L] <= 0.5) # the x1 > 0.5 half-space carries no weight
kept <- w > 0
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  updateState = FALSE,
  seed = 11L
)

# 300 forests of 20 trees each, all from the C-side stream the control seed
# fixes, so the counts are deterministic
drawAndCountLeaves <- function(weights) {
  sampler <- dbarts::dbarts(x, y, weights = weights, control = control)
  empty <- 0L
  weightless <- 0L
  leaves <- 0L
  structure <- character(300L)
  for (i in seq_len(300L)) {
    sampler$sampleTreesFromPrior()
    members <- sampler$getTrees(current = TRUE)
    weighted <- sampler$getTrees(
      current = TRUE,
      newdata = x[kept, , drop = FALSE]
    )
    isLeaf <- members$var == -1L
    leaves <- leaves + sum(isLeaf)
    empty <- empty + sum(members$n[isLeaf] == 0L)
    weightless <- weightless + sum(weighted$n[isLeaf] == 0L)
    structure[i] <- paste0(members$n, ":", members$var, collapse = "|")
  }
  list(
    leaves = leaves,
    empty = empty,
    weightless = weightless,
    structure = structure
  )
}

weighted <- drawAndCountLeaves(w)
expect_true(weighted$leaves > 6000L) # the trees grew; the check bites
expect_equal(weighted$empty, 0L)
# leaves that no positive-weight row reaches are legal, and common
expect_true(weighted$weightless > 100L)

# the prior over trees does not depend on the weights
counted <- drawAndCountLeaves(NULL)
expect_identical(weighted$structure, counted$structure)

# nor on an active-row mask, an all-zeros one included: the forests are the
# same draws, not one bare root per tree
masked <- dbarts::dbarts(x, y, control = control)
masked$setActiveRows(rep(0, n))
masked$sampleTreesFromPrior()
maskedTrees <- masked$getTrees(current = TRUE)
expect_identical(
  paste0(maskedTrees$n, ":", maskedTrees$var, collapse = "|"),
  counted$structure[1L]
)

rm(
  n,
  x,
  w,
  kept,
  y,
  control,
  drawAndCountLeaves,
  weighted,
  counted,
  masked,
  maskedTrees
)
