# sampleTreesFromPrior draws the tree prior CONDITIONED on the empty-leaf-free
# set the move kernels price, by per-tree
# rejection. The conditioning predicate is the veto's own - positive WEIGHT,
# not membership - so under a zero-weight half-space no leaf of a from-prior
# forest may hold only zero-weight rows. The projection this replaced collapsed
# on the member count, which spares such a leaf; a birth out of one then
# compared -HUGE_VAL to -HUGE_VAL, a NaN ratio that rejected silently.
#
# The oracle needs no tree walk: routing ONLY the positive-weight rows through
# the drawn trees (getTrees(newdata = )) counts, per node, the rows the
# likelihood can see, so the law says every leaf reports n > 0.

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
# fixes, so the count is deterministic
drawAndCountEmptyLeaves <- function(weights) {
  args <- list(x, y, control = control)
  if (!is.null(weights)) {
    args$weights <- weights
  }
  sampler <- suppressWarnings(do.call(dbarts::dbarts, args))
  empty <- 0L
  leaves <- 0L
  for (i in seq_len(300L)) {
    sampler$sampleTreesFromPrior()
    trees <- sampler$getTrees(
      current = TRUE,
      newdata = x[kept, , drop = FALSE]
    )
    bottoms <- trees[trees$var == -1L, ]
    leaves <- leaves + nrow(bottoms)
    empty <- empty + sum(bottoms$n == 0L)
  }
  c(leaves = leaves, empty = empty)
}

weighted <- drawAndCountEmptyLeaves(w)
expect_true(weighted[["leaves"]] > 6000L) # the trees grew; the check bites
expect_equal(weighted[["empty"]], 0L)

# non-vacuity: with no weight vector installed the same fixture conditions on
# membership alone, and leaves that no positive-weight row reaches are common
counted <- drawAndCountEmptyLeaves(NULL)
expect_true(counted[["empty"]] > 100L)

rm(n, x, w, kept, y, control, drawAndCountEmptyLeaves, weighted, counted)
