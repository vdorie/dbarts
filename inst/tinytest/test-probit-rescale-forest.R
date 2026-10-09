# The probit rescaling step, from R: that it sits ahead of every recorded
# channel, runs per chain on any number of threads, runs finite under an
# offset, a mask and an infinite prior scale, and that a sampler's state
# carries what it needs. The step's conditional, its mapping, that it acts on
# every chain and its declines are tests/cpp gates; its posterior is the
# probit-k-scale-exact gate; which fits take it is the equivalence baseline's
# record.

set.seed(99)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
eta <- 1.5 * sin(3 * x[, 1L]) + x[, 2L] - 0.7
yBinary <- as.numeric(eta + rnorm(n) > 0)
xTest <- matrix(runif(30L), 10L, 3L, dimnames = list(NULL, colnames(x)))

stepControl <- function(n.chains = 1L, n.threads = 1L) {
  dbarts::dbartsControl(
    n.trees = 10L,
    n.chains = n.chains,
    n.threads = n.threads,
    n.burn = 0L,
    n.samples = 20L,
    seed = 7L
  )
}

# ---- per chain, on any number of threads ----

twoChains <- function(n.threads) {
  dbarts::dbarts(
    x,
    yBinary,
    control = stepControl(n.chains = 2L, n.threads = n.threads)
  )$run(10L, 20L)
}
oneThread <- twoChains(1L)
expect_identical(oneThread, twoChains(2L))
rm(oneThread)

# ---- ahead of every recorded channel ----

# the saved trees replayed reproduce each kept draw's training and test fits,
# which holds only if the step rescales the leaves before any of them is
# written
treesFit <- dbarts::bart(
  x,
  yBinary,
  test = xTest,
  n.trees = 10L,
  n.samples = 30L,
  n.burn = 30L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  seed = 7L,
  verbose = FALSE
)
expect_true(
  max(abs(predict(treesFit, x, type = "bart") - treesFit$yhat.train)) < 1e-10
)
expect_true(
  max(abs(predict(treesFit, xTest, type = "bart") - treesFit$yhat.test)) < 1e-10
)
rm(treesFit)

# ---- an offset, a mask installed mid-run and an infinite prior scale ----

masked <- dbarts::dbarts(
  x,
  yBinary,
  offset = 0.3 * x[, 3L],
  leaf.prior = normal(k = chi(1.5, Inf)),
  control = stepControl()
)
invisible(masked$run(20L, 1L))
masked$setActiveRows(rep(c(1, 1, 1, 1, 0), length.out = n))
maskedRun <- masked$run(20L, 20L)
expect_true(all(is.finite(maskedRun$train)) && all(is.finite(maskedRun$k)))
rm(masked, maskedRun)

# ---- the state round-trip ----

# a sampler restored from the store and a reload of the same store continue
# identically, k and the fits both
stored <- dbarts::dbarts(x, yBinary, control = stepControl())
invisible(stored$run(10L, 1L))
stored$storeState()
serialized <- tempfile(fileext = ".rds")
saveRDS(stored, serialized)
restored <- dbarts::dbarts(x, yBinary, control = stepControl())
restored$setState(stored$state)
reloaded <- readRDS(serialized)
unlink(serialized)
expect_identical(restored$run(0L, 10L), reloaded$run(0L, 10L))
rm(stored, restored, reloaded, serialized)
