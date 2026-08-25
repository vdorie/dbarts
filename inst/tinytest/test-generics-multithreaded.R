# predict's n.threads across the whole predict surface. The engine partitions
# the saved-tree replay by (chain, draw), each partition writing its own rows
# and reducing nothing across workers, so the ONLY thing a thread count may
# change is the fan-out. Both halves are pinned here: that the answer is
# identical bit for bit at every count (identical(), not expect_equal), and
# that the count actually reached the engine and dealt the slabs out - read off
# a deterministic worker/partition channel rather than a timing ratio, which
# would be a benchmark rather than a test.
#
# Every fit and every predict here bounds itself at 2 threads.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
friedman <- testData
rm(testData)

# A replay this small is far below the traversal cutoff, so it would run inline
# whatever count it was given and every identity below would be vacuous.
# Override the cutoff for this file and restore it at the end.
previousCutoff <- .Call(dbarts:::C_dbarts_bartcore_setPredictParallelCutoff, 1L)

# The channel: the count the engine resolved, the workers it ran, and the
# slab -> worker map, all recorded before any worker started. Reported as a
# named logical rather than as expectations, so a failure names WHICH property
# broke while the expectation itself stays at the top level, where tinytest
# registers it.
partition <- function(n.threads, n.slabs) {
  part <- .Call(dbarts:::C_dbarts_bartcore_lastPredictPartition)
  workers <- min(n.threads, n.slabs)
  c(
    resolved = identical(part$resolved, n.threads),
    workers = identical(part$n.workers, workers),
    slabs = identical(length(part$worker), n.slabs),
    # a contiguous block partition covering every slab exactly once: every
    # worker owns at least one slab, the map never goes backwards, and the
    # block sizes differ by at most one
    labelled = identical(sort(unique(part$worker)), seq_len(workers) - 1L),
    blocked = !is.unsorted(part$worker),
    balanced = max(table(part$worker)) - min(table(part$worker)) <= 1L
  )
}
partitionOK <- c(
  resolved = TRUE,
  workers = TRUE,
  slabs = TRUE,
  labelled = TRUE,
  blocked = TRUE,
  balanced = TRUE
)

n.samples <- 5L

# ---- gaussian, one chain: 5 slabs over 2 workers, so the block sizes are 3
# and 2 and an even-division assumption in the partition would drop a slab ----
gaussianFit <- dbarts::bart(
  friedman$x,
  friedman$y,
  ndpost = n.samples,
  nskip = 5L,
  ntree = 5L,
  nthread = 1L,
  verbose = FALSE,
  keeptrees = TRUE
)
serial <- predict(gaussianFit, friedman$x, n.threads = 1L)
expect_identical(partition(1L, n.samples), partitionOK)
threaded <- predict(gaussianFit, friedman$x, n.threads = 2L)
expect_identical(partition(2L, n.samples), partitionOK)
expect_identical(serial, threaded)
# the replay is still the recorded fit, thread count or not
expect_equal(threaded, gaussianFit$yhat.train)

# ---- gaussian, two chains: 10 slabs, and the chain axis sits inside the
# partition rather than being the partition ----
chainFit <- dbarts::bart(
  friedman$x,
  friedman$y,
  ndpost = n.samples,
  nskip = 5L,
  ntree = 5L,
  nchain = 2L,
  nthread = 1L,
  verbose = FALSE,
  keeptrees = TRUE
)
serialChains <- predict(chainFit, friedman$x, n.threads = 1L)
expect_identical(serialChains, predict(chainFit, friedman$x, n.threads = 2L))
expect_identical(partition(2L, 2L * n.samples), partitionOK)
expect_equal(serialChains, chainFit$yhat.train)

# a sparse test source routes its rows resident and replays over the same
# partition
if (requireNamespace("Matrix", quietly = TRUE)) {
  sparse.x <- Matrix::Matrix(friedman$x, sparse = TRUE)
  expect_identical(
    predict(gaussianFit, sparse.x, n.threads = 1L),
    predict(gaussianFit, sparse.x, n.threads = 2L)
  )
  expect_identical(partition(2L, n.samples), partitionOK)
  rm(sparse.x)
}

# ---- binary: the latent arm and the probability arm both replay ----
set.seed(6001)
nBin <- 120L
xBin <- matrix(rnorm(nBin * 3L), nBin, 3L)
yBin <- rbinom(nBin, 1L, pnorm(xBin %*% c(0.12, -0.05, 0.3)))
binaryFit <- dbarts::bart2(
  xBin,
  yBin,
  n.trees = 5L,
  n.samples = n.samples,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
for (predictType in c("bart", "ev")) {
  expect_identical(
    predict(binaryFit, xBin, type = predictType, n.threads = 1L),
    predict(binaryFit, xBin, type = predictType, n.threads = 2L)
  )
  expect_identical(partition(2L, n.samples), partitionOK)
}

# ---- multinomial: K forests per slab, softmaxed inside the slab, so the
# partition must keep a (slab, row) pair whole ----
set.seed(2201)
nMulti <- 120L
xMulti <- matrix(runif(nMulti * 3L), nMulti, 3L)
etaMulti <- cbind(
  2 * (xMulti[, 1L] - 0.5),
  xMulti[, 2L] - xMulti[, 3L],
  1.5 * (xMulti[, 3L] - 0.5)
)
probsMulti <- exp(etaMulti) / rowSums(exp(etaMulti))
yMulti <- factor(
  c("lo", "mid", "hi")[vapply(
    seq_len(nMulti),
    function(i) sample.int(3L, 1L, prob = probsMulti[i, ]),
    integer(1L)
  )],
  levels = c("lo", "mid", "hi")
)
multinomialFit <- dbarts::bart2(
  xMulti,
  yMulti,
  family = "multinomial",
  n.trees = 5L,
  n.samples = n.samples,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(
  predict(multinomialFit, xMulti, n.threads = 1L),
  predict(multinomialFit, xMulti, n.threads = 2L)
)
expect_identical(partition(2L, n.samples), partitionOK)

# ---- heteroscedastic: the mean and the s(x) surface are two replays, each
# with its own partition, returned together ----
set.seed(717)
nVar <- 150L
xVar <- matrix(runif(nVar), nVar, 1L)
yVar <- 2 * xVar[, 1L] + ifelse(xVar[, 1L] < 0.5, 0.3, 1.5) * rnorm(nVar)
varianceFit <- dbarts::bart2(
  xVar,
  yVar,
  variance = dbarts::varianceForest(n.trees = 5L),
  n.trees = 5L,
  n.samples = n.samples,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
serialVar <- predict(varianceFit, xVar, type = "bart", n.threads = 1L)
threadedVar <- predict(varianceFit, xVar, type = "bart", n.threads = 2L)
expect_identical(serialVar, threadedVar)
# s(x) rides back as an attribute, and it is the SECOND replay: an identity on
# the mean alone would never reach the variance entry point
expect_false(is.null(attr(serialVar, "s")))
expect_identical(attr(serialVar, "s"), attr(threadedVar, "s"))
expect_identical(partition(2L, n.samples), partitionOK)

# ---- amplitude-coupled (BCF glue): the per-forest replay AND the plain
# predict that routes through it. A blend reaches the sampler only here, so
# without this leg the coupling's predict would stay untested ----
set.seed(811)
nBcf <- 60L
xBcf <- matrix(runif(nBcf * 3L), nBcf, 3L)
colnames(xBcf) <- paste0("x", seq_len(3L))
zBcf <- rbinom(nBcf, 1L, 0.5)
yBcf <- 2 *
  sin(pi * xBcf[, 1L]) +
  zBcf * (1 + 2 * xBcf[, 2L]) +
  rnorm(nBcf, 0, 0.2)

bcfSampler <- dbarts::dbarts(
  xBcf,
  yBcf,
  forests = list(dbarts::forest(), dbarts::forest(basis = ~ factor(zBcf))),
  control = dbarts::dbartsControl(
    n.threads = 1L,
    n.trees = 8L,
    n.chains = 1L,
    n.samples = n.samples,
    n.burn = 2L,
    keepTrees = FALSE,
    updateState = FALSE,
    seed = 811L
  )
)
bcfBurn <- dbarts:::runWithBurnIn(bcfSampler, bcfSampler$control, TRUE)
bcfFit <- dbarts:::packageBartResults(
  bcfSampler,
  bcfBurn$samples,
  bcfBurn$burnInSigma,
  bcfBurn$burnInK,
  TRUE,
  TRUE
)

expect_identical(
  predict(bcfFit, xBcf, type = "forest", n.threads = 1L),
  predict(bcfFit, xBcf, type = "forest", n.threads = 2L)
)
expect_identical(partition(2L, n.samples), partitionOK)
expect_identical(
  predict(bcfFit, xBcf, type = "bart", bases = bcfFit$bases, n.threads = 1L),
  predict(bcfFit, xBcf, type = "bart", bases = bcfFit$bases, n.threads = 2L)
)
expect_identical(partition(2L, n.samples), partitionOK)

# ---- rbart_vi ----
set.seed(919)
nRe <- 80L
xRe <- matrix(runif(nRe * 2L), nRe, 2L)
g <- factor(rep_len(seq_len(4L), nRe))
yRe <- xRe[, 1L] + rnorm(4L)[g] + rnorm(nRe, 0, 0.3)
rbartFit <- dbarts::rbart_vi(
  xRe,
  yRe,
  group.by = g,
  n.trees = 5L,
  n.samples = n.samples,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(
  predict(rbartFit, xRe, g, n.threads = 1L),
  predict(rbartFit, xRe, g, n.threads = 2L)
)
# rbart's Gibbs loop advances its sampler one draw at a time, so its store
# holds a single slab: the worker count is clamped to it rather than spawning
# an idle second worker
expect_identical(partition(2L, 1L), partitionOK)
# the formal is appended last, so a caller's THIRD positional argument is
# still group.by
expect_identical(
  predict(rbartFit, xRe, g, n.threads = 2L),
  predict(rbartFit, xRe, group.by = g, n.threads = 2L)
)

# ---- the formal is LAST on every predict method, immediately before '...' ----
for (method in c(
  "predict.bart",
  "predict.bartMultinomial",
  "predict.bartOrdinal",
  "predict.bartNegbin",
  "predict.bartHurdle",
  "predict.rbart"
)) {
  formalNames <- names(formals(utils::getFromNamespace(method, "dbarts")))
  expect_identical(
    formalNames[c(length(formalNames) - 1L, length(formalNames))],
    c("n.threads", "...")
  )
}

# ---- refusals, by name and echoing the value ----
for (bad in list(0L, -1L, NA_integer_, "two", c(1L, 2L), TRUE, NULL)) {
  expect_error(
    predict(gaussianFit, friedman$x, n.threads = bad),
    "'n.threads' must be a single positive integer"
  )
}
expect_error(predict(gaussianFit, friedman$x, n.threads = 0L), "not 0L")
expect_error(predict(multinomialFit, xMulti, n.threads = -3L), "not -3L")
# the refusal precedes the type = "forest" and blend returns, so an
# amplitude-coupled fit's value is checked too rather than being carried into
# an arm that validates nothing
expect_error(
  predict(bcfFit, xBcf, type = "forest", n.threads = 0L),
  "'n.threads' must be a single positive integer"
)
expect_error(
  predict(bcfFit, xBcf, type = "bart", bases = bcfFit$bases, n.threads = NA),
  "'n.threads' must be a single positive integer"
)

# ---- the bridge refuses a non-positive count on its own, so a caller that
# reaches the entry point directly is never left with an unwritten buffer ----
expect_error(
  .Call(
    dbarts:::C_dbarts_bartcore_predict,
    gaussianFit$fit$getPointer(),
    friedman$x,
    NULL,
    0L
  ),
  "number of threads"
)

invisible(.Call(
  dbarts:::C_dbarts_bartcore_setPredictParallelCutoff,
  previousCutoff
))

rm(
  friedman,
  previousCutoff,
  partition,
  partitionOK,
  n.samples,
  gaussianFit,
  serial,
  threaded,
  chainFit,
  serialChains,
  nBin,
  xBin,
  yBin,
  binaryFit,
  predictType,
  nMulti,
  xMulti,
  etaMulti,
  probsMulti,
  yMulti,
  multinomialFit,
  nVar,
  xVar,
  yVar,
  varianceFit,
  serialVar,
  threadedVar,
  nBcf,
  xBcf,
  zBcf,
  yBcf,
  bcfSampler,
  bcfBurn,
  bcfFit,
  nRe,
  xRe,
  yRe,
  g,
  rbartFit,
  method,
  formalNames,
  bad
)
