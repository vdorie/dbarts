# pdbart/pd2bart obtain their sampler through bart(sampleronly = TRUE), which
# hands one back unrun, and then run only the burn-in phase. Their keepTrees
# branches predict from the saved-tree store instead of running, so without a
# sampling phase they read a store no draw was ever recorded into: the
# partial-dependence surface came back flat (sd 0 at every level), silently.
# The sampling phase now happens where the burn-in does, so the two branches
# compute the same surface from the same draws.

set.seed(11L)
n <- 80L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("x1", "x2", "x3")
y <- 3 * x[, 1L] - 2 * x[, 2L] + rnorm(n, 0, 0.4)
levels.pd <- list(c(0.2, 0.5, 0.8))

pdFrom <- function(keeptrees) {
  set.seed(5L)
  dbarts::pdbart(
    x,
    y,
    xind = 1L,
    levs = levels.pd,
    keeptrees = keeptrees,
    nskip = 0L,
    ndpost = 6L,
    ntree = 10L,
    pl = FALSE,
    verbose = FALSE
  )
}

pd.kept <- pdFrom(TRUE)
pd.run <- pdFrom(FALSE)

# the surface varies over the draws and over the levels: the defect returned a
# single value everywhere
expect_equal(dim(pd.kept$fd[[1L]]), c(6L, 3L))
expect_true(all(apply(pd.kept$fd[[1L]], 2L, sd) > 1e-8))
expect_true(diff(range(colMeans(pd.kept$fd[[1L]]))) > 1e-8)

# and it is the same surface the non-keepTrees branch computes by re-running:
# same seed, same draws, and a saved-tree replay reproduces the run's own test
# fits exactly (measured 0 here, held to 1e-12 across platforms)
expect_true(max(abs(pd.kept$fd[[1L]] - pd.run$fd[[1L]])) < 1e-12)

# pd2bart shares the prologue and so is fixed with it
set.seed(5L)
pd2.kept <- dbarts::pd2bart(
  x,
  y,
  xind = 1:2,
  levs = list(c(0.2, 0.8), c(0.2, 0.8)),
  keeptrees = TRUE,
  nskip = 0L,
  ndpost = 6L,
  ntree = 10L,
  pl = FALSE,
  verbose = FALSE
)
expect_equal(dim(pd2.kept$fd), c(6L, 4L))
expect_true(all(apply(pd2.kept$fd, 2L, sd) > 1e-8))

rm(n, x, y, levels.pd, pdFrom, pd.kept, pd.run, pd2.kept)
