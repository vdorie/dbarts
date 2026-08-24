source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

makeSampler <- function(n.trees = 25L) {
  dbarts::dbarts(
    x,
    y,
    control = dbarts::dbartsControl(
      n.trees = n.trees,
      n.chains = 1L,
      updateState = FALSE,
      keepTrees = FALSE
    )
  )
}

## grow-from-root in place yields a well-formed forest with finite predictions
sampler <- makeSampler()
sampler$growFromRoot(3L)
grownFit <- as.vector(sampler$predict(x))
expect_true(all(is.finite(grownFit)))

## the exact sampler continues from the grown forest with finite fits
res <- sampler$run(0L, 10L)
expect_true(all(is.finite(res$train)))

## seeded reproducibility: same seed -> identical grown fit
set.seed(42L)
s1 <- makeSampler()
s1$growFromRoot(2L)
f1 <- as.vector(s1$predict(x))
set.seed(42L)
s2 <- makeSampler()
s2$growFromRoot(2L)
f2 <- as.vector(s2$predict(x))
expect_equal(f1, f2)

## a different seed grows a different forest (grow consumes RNG)
set.seed(7L)
s3 <- makeSampler()
s3$growFromRoot(2L)
f3 <- as.vector(s3$predict(x))
expect_false(isTRUE(all.equal(f1, f3)))

## cross-sampler workflow: donor$growFromRoot -> target$installTrees round trip
## reproduces the grown forest with no bespoke install code
donor <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
donor$growFromRoot(2L)
target <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
target$installTrees(donor)
expect_equal(as.vector(target$predict(x)), as.vector(donor$predict(x)))

## n.sweeps must be a single positive integer
expect_error(makeSampler()$growFromRoot(0L), "positive integer")
expect_error(makeSampler()$growFromRoot(-2L), "positive integer")

## the linear-leaf model refuses grow-from-root with an informative error
linSampler <- dbarts::dbarts(
  x,
  y,
  node.prior = linear(1L),
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE
  )
)
expect_error(linSampler$growFromRoot(2L), "constant-leaf")

## updateState = NA (the default) respects control@updateState, matching the
## sibling initializers sampleTreesFromPrior/sampleNodeParametersFromPrior. An
## updateState = FALSE control stores nothing on the default, and only an
## explicit updateState = TRUE opts in; an updateState = TRUE control's default
## refreshes the cached state (the pre-fix default of FALSE would have skipped
## that store, leaving the realized cache stale).
source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

nostore <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 10L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
nostore$growFromRoot(2L) # default updateState = NA, control FALSE -> no store
expect_null(nostore$state)

optin <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 10L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
optin$growFromRoot(2L, updateState = TRUE) # explicit TRUE stores anyway
expect_false(is.null(optin$state))

store <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 10L,
    n.chains = 1L,
    updateState = TRUE,
    keepTrees = FALSE
  )
)
initialState <- store$state # realize the cache on the initial (stump) forest
store$growFromRoot(2L) # default NA, control TRUE -> stores, refreshing state
expect_false(is.null(store$state))
expect_false(statesAgree(store$state, initialState, expect = FALSE))

rm(nostore, optin, store, initialState)

## a heteroscedastic sampler grows against the variance-weighted precisions.
## Under a variance forest the global sigma is pinned at 1 on the working
## scale, so a scan reading unit precisions prices every split against a
## residual variance the data does not have. The homoscedastic control is
## the reference: its sigma carries the scale, so its own scan was always
## right, and it is fit to the same data by the same recipe. Reweighting by
## 1/s^2(x) is more efficient
## where the noise varies 900-fold in variance, so the heteroscedastic init
## must beat it - and does, once the pre-step is there (cor 0.96 vs the
## control's 0.86; without the pre-step it was 0.14, and the quiet region's
## RMSE 1.32 against the control's 0.55).
set.seed(53, sample.kind = "Rejection")
nGrow <- 400L
xGrow <- cbind(runif(nGrow), runif(nGrow))
muGrow <- 4 * xGrow[, 1L]
yGrow <- muGrow + ifelse(xGrow[, 2L] < 0.5, 0.1, 3) * rnorm(nGrow)
growControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 25L,
  n.burn = 0L,
  verbose = FALSE
)
growInit <- function(sampler) {
  invisible(sampler$run(50L, 0L)) # fit the scale surface the scan then reads
  sampler$growFromRoot(1L)
  as.vector(sampler$predict(xGrow))
}
fitHeteroscedastic <- growInit(dbarts::dbarts(
  xGrow,
  yGrow,
  control = growControl,
  variance = dbarts::varianceForest(n.trees = 10L)
))
fitHomoscedastic <- growInit(dbarts::dbarts(
  xGrow,
  yGrow,
  control = growControl
))
expect_true(
  cor(fitHeteroscedastic, muGrow) > cor(fitHomoscedastic, muGrow)
)
## and the gain is where the weighting puts it: the low-noise half, which the
## precisions promote from a 900th of the weight to parity
quiet <- xGrow[, 2L] < 0.5
rmseQuiet <- function(fit) sqrt(mean((fit - muGrow)[quiet]^2))
expect_true(rmseQuiet(fitHeteroscedastic) < 0.5 * rmseQuiet(fitHomoscedastic))

rm(
  nGrow,
  xGrow,
  muGrow,
  yGrow,
  growControl,
  growInit,
  fitHeteroscedastic,
  fitHomoscedastic,
  quiet,
  rmseQuiet
)

## Missing ordinal predictors: the grow scan enumerates BOTH missing directions
## of every cut where a node holds missing rows and scores each with those rows
## in the child it sends them to, so a rule the initializer places on such a
## column carries a real routing decision and every row still lands in a leaf.
set.seed(20260824L)
nMissing <- 200L
xMissing <- matrix(runif(nMissing * 3L), nMissing, 3L)
missingRows <- seq.int(1L, nMissing, by = 5L)
## missingness that carries the signal: the routing decision has to be worth
## making for the scan to be exercised where it matters
missingSignal <- ifelse(seq_len(nMissing) %in% missingRows, 4, 0)
yMissing <- 3 * xMissing[, 1L] + missingSignal + rnorm(nMissing, 0, 0.2)
xMissing[missingRows, 1L] <- NA_real_

missingSampler <- dbarts::dbarts(
  xMissing,
  yMissing,
  control = dbarts::dbartsControl(
    n.trees = 25L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
missingSampler$growFromRoot(2L)
missingTrees <- missingSampler$getTrees()

## the initializer split on the NA-bearing column, and every such rule names a
## direction for its missing rows
onColumn1 <- missingTrees$var == 1L
expect_true(any(onColumn1))
expect_true(all(missingTrees$missing[onColumn1] %in% c("L", "R")))
## and the grown forest is a legal chain state on data the routing reaches:
## no empty leaf, every row predicted
expect_true(all(missingTrees$n[missingTrees$var == -1L] > 0L))
expect_true(all(is.finite(missingSampler$predict(xMissing))))

rm(
  nMissing,
  xMissing,
  missingRows,
  missingSignal,
  yMissing,
  missingSampler,
  missingTrees,
  onColumn1
)
