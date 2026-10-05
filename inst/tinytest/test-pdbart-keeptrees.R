# pdbart and pd2bart predict from the saved trees of the fit they make, so the
# sampling phase has to record its draws even with no burn-in: a run that
# recorded none returned a flat surface (sd 0 at every level), silently. A
# BayesTree-spelled call with keeptrees = TRUE and nskip = 0 is translated
# and fits the same model as the call in bart's names.

onceState <- dbarts:::onceWarnState
savedKeys <- grep("^tombstone\\.pd2?bart", ls(onceState), value = TRUE)
savedState <- mget(savedKeys, envir = onceState)

set.seed(11L)
n <- 80L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("x1", "x2", "x3")
y <- 3 * x[, 1L] - 2 * x[, 2L] + rnorm(n, 0, 0.4)
levels.pd <- list(c(0.2, 0.5, 0.8))

pd.translated <- suppressWarnings(dbarts::pdbart(
  x,
  y,
  xind = 1L,
  levs = levels.pd,
  keeptrees = TRUE,
  nskip = 0L,
  ndpost = 6L,
  ntree = 10L,
  nchain = 1L,
  nthread = 1L,
  seed = 5L,
  pl = FALSE,
  verbose = FALSE
))
pd.named <- dbarts::pdbart(
  x,
  y,
  xind = 1L,
  levs = levels.pd,
  keepTrees = TRUE,
  n.burn = 0L,
  n.samples = 6L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  pl = FALSE,
  verbose = FALSE
)

# the surface varies over the draws and over the levels: the defect returned a
# single value everywhere
expect_equal(dim(pd.translated$fd[[1L]]), c(6L, 3L))
expect_true(all(apply(pd.translated$fd[[1L]], 2L, sd) > 1e-8))
expect_true(diff(range(colMeans(pd.translated$fd[[1L]]))) > 1e-8)
expect_identical(pd.translated$fd, pd.named$fd)

# pd2bart shares the prologue and so is fixed with it
pd2.kept <- dbarts::pd2bart(
  x,
  y,
  xind = 1:2,
  levs = list(c(0.2, 0.8), c(0.2, 0.8)),
  n.burn = 0L,
  n.samples = 6L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 5L,
  pl = FALSE,
  verbose = FALSE
)
expect_equal(dim(pd2.kept$fd), c(6L, 4L))
expect_true(all(apply(pd2.kept$fd, 2L, sd) > 1e-8))

for (key in grep("^tombstone\\.pd2?bart", ls(onceState), value = TRUE)) {
  onceState[[key]] <- if (key %in% savedKeys) savedState[[key]]
}
rm(n, x, y, levels.pd, pd.translated, pd.named, pd2.kept)
rm(key, onceState, savedKeys, savedState)
