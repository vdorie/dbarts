source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

df <- with(testData, data.frame(x, y))

# test that base bart extracts trees correctly
n.trees <- 3L
n.samples <- 4L
fit <- dbarts::bart(
  y ~ .,
  df,
  nthread = 1L,
  ntree = n.trees,
  nskip = 0L,
  ndpost = n.samples,
  keeptrees = TRUE,
  verbose = FALSE
)
allTrees <- dbarts::extract(fit, "trees")

# a single-forest fit has no forest column, and extract always carries chain
expect_equal(
  colnames(allTrees),
  c("chain", "sample", "tree", "n", "var", "value")
)
expect_true(all(allTrees$chain == 1L))

combinations <- data.frame(
  sample = rep(seq_len(n.samples), each = n.trees),
  tree = rep(seq_len(n.trees), times = n.samples)
)
expect_true(
  all(
    paste0(combinations$sample, ";", combinations$tree) %in%
      paste0(allTrees$sample, ";", allTrees$tree)
  )
)

individualSamples <- lapply(
  seq_len(n.samples),
  function(i) extract(fit, "trees", sampleNums = i)
)
individualSamples <- Reduce(rbind, individualSamples)
row.names(individualSamples) <- as.character(seq_len(nrow(individualSamples)))

expect_equal(allTrees, individualSamples)

# extract's own formals (sample, combineChains, contribution) do not reach
# getTrees (bart.Rd's 'Extracting Trees' section documents chainNums/
# sampleNums/treeNums/newdata/forest as accepted there); each is refused by
# name instead of silently corrupting the call or partial-matching one of
# getTrees's differently-named formals (sample -> sampleNums). 'forest' is
# getTrees' own formal name and forwards instead of being refused.
treesArgReason <- function(arg) {
  paste0(
    "'",
    arg,
    "' is not used when type = \"trees\"; the sampler's getTrees ",
    "accepts 'chainNums', 'sampleNums', 'treeNums', 'current', 'newdata', ",
    "and 'forest' instead (see 'Extracting Trees' in ?bart)"
  )
}
expect_error(
  extract(fit, type = "trees", sample = 1L),
  treesArgReason("sample"),
  fixed = TRUE
)
expect_error(
  extract(fit, type = "trees", sample = "train"),
  treesArgReason("sample"),
  fixed = TRUE
)
expect_error(
  extract(fit, "trees", "train"),
  treesArgReason("sample"),
  fixed = TRUE
)
expect_error(
  extract(fit, type = "trees", combineChains = FALSE),
  treesArgReason("combineChains"),
  fixed = TRUE
)
expect_error(
  extract(fit, type = "trees", contribution = TRUE),
  treesArgReason("contribution"),
  fixed = TRUE
)

# 'forest' forwards to getTrees rather than being refused: a single-forest
# fit's forest = 1L read is bitwise its default, and an out-of-range or
# malformed forest is refused exactly as $getTrees refuses it directly
expect_identical(
  extract(fit, type = "trees", forest = 1L),
  allTrees
)
expect_error(
  extract(fit, type = "trees", forest = 2L),
  "out of range"
)
expect_error(
  extract(fit, type = "trees", forest = 0L),
  "'forest' must be a single positive integer"
)

rm(individualSamples, combinations, allTrees, fit, n.samples, n.trees)

## ---------------------------------------------------------------------------
## extract(type = "trees") on a keepSampler-only fit (keeptrees FALSE,
## keepSampler TRUE): follows plotTree's own documented fallback, the
## sampler's CURRENT trees - no sample column - rather than a saved
## history, since keeptrees FALSE means no history was kept.
## ---------------------------------------------------------------------------
n.trees <- 3L
fitKeepSampler <- dbarts::bart(
  y ~ .,
  df,
  nthread = 1L,
  ntree = n.trees,
  nskip = 0L,
  ndpost = 4L,
  keeptrees = FALSE,
  keepsampler = TRUE,
  verbose = FALSE
)
currentTrees <- dbarts::extract(fitKeepSampler, "trees")
expect_equal(colnames(currentTrees), c("chain", "tree", "n", "var", "value"))
expect_true(nrow(currentTrees) > 0L)
expect_equal(
  currentTrees[colnames(currentTrees) != "chain"],
  fitKeepSampler$fit$getTrees()
)

fitKept <- dbarts::bart(
  y ~ .,
  df,
  nthread = 1L,
  ntree = n.trees,
  nskip = 0L,
  ndpost = 4L,
  keeptrees = TRUE,
  verbose = FALSE
)
keptTrees <- dbarts::extract(fitKept, "trees")
expect_equal(
  colnames(keptTrees),
  c("chain", "sample", "tree", "n", "var", "value")
)
expect_true(nrow(keptTrees) > 0L)

# an EMPTY index selection is a zero-row answer rather than an error, through
# both entrances: the bridge emits each column from a gather that holds no
# block at all
expect_equal(
  nrow(dbarts::extract(fitKept, "trees", treeNums = integer(0))),
  0L
)
expect_equal(
  nrow(fitKeepSampler$fit$getTrees(treeNums = integer(0))),
  0L
)

rm(n.trees, fitKeepSampler, currentTrees, fitKept, keptTrees)


## ---------------------------------------------------------------------------
## getTrees(current = TRUE) returns the live working trees even for a keepTrees
## sampler: no sample dimension, identical to what a keepTrees = FALSE sampler
## reports, and a valid live oracle after a partial update (unlike the saved
## snapshots, whose n replays the current predictor through frozen structure).
## ---------------------------------------------------------------------------
set.seed(7L)
n <- 50L
x <- rnorm(n)
y <- x + rnorm(n)
makeSamplerTrees <- function(keepTrees) {
  ctrl <- dbarts::dbartsControl(
    n.chains = 1L,
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = TRUE,
    keepTrees = keepTrees,
    verbose = FALSE,
    seed = 9L
  )
  sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
  invisible(sampler$run(30L, 5L))
  sampler
}
keptSampler <- makeSamplerTrees(TRUE)
liveSampler <- makeSamplerTrees(FALSE)

saved <- keptSampler$getTrees() # snapshots
current <- keptSampler$getTrees(current = TRUE) # live working trees

expect_true("sample" %in% names(saved)) # snapshots carry a sample dimension
expect_false("sample" %in% names(current)) # live trees do not
# the live trees are exactly what an otherwise-identical keepTrees = FALSE
# sampler reports
expect_equal(current, liveSampler$getTrees())
# current = TRUE is ignored (harmlessly) when there are no saved trees
expect_equal(liveSampler$getTrees(current = TRUE), liveSampler$getTrees())
# and current = TRUE takes the sample handling of a keepTrees = FALSE sampler:
# the live trees carry no sample dimension, so a sampleNums filter names
# nothing to filter and is dropped with a warning rather than applied to the
# saved store or range-checked against it
warnings.currentSample <- captureWarnings(
  currentFiltered <- keptSampler$getTrees(current = TRUE, sampleNums = 1L)
)
expect_equal(length(warnings.currentSample), 1L)
expect_true(any(grepl(
  "sampleNums ignored if current is TRUE",
  vapply(warnings.currentSample, conditionMessage, ""),
  fixed = TRUE
)))
expect_true(any(vapply(
  warnings.currentSample,
  inherits,
  logical(1L),
  "dbartsIgnoredArgWarning"
)))
expect_equal(currentFiltered, current)

# a partial update keeps the live trees valid; getTrees(current = TRUE) sees it
inst <- keptSampler$setPredictor(
  x + rnorm(n, sd = 0.5),
  "x",
  forceUpdate = "partial"
)
liveTrees <- keptSampler$getTrees(current = TRUE)
expect_false(any(liveTrees$var == -1L & liveTrees$n == 0L))
expect_true(length(inst) == n)

rm(keptSampler, liveSampler, makeSamplerTrees, saved, current, liveTrees, inst)
rm(warnings.currentSample, currentFiltered)
rm(x, y, n)


## ---------------------------------------------------------------------------
## getTrees(newdata = X) routes X through each tree so 'n' counts that data.
## Routing the training design matrix reproduces the default counts on every
## path (saved snapshots, live-from-keepTrees, live), and an independent R-side
## descent reproduces the counts for arbitrary data.
## ---------------------------------------------------------------------------
set.seed(11L)
n <- 60L
x <- rnorm(n)
y <- x + rnorm(n)
ctrl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 0L,
  updateState = TRUE,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 3L
)
keptSampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(keptSampler$run(30L, 5L))

ctrl@keepTrees <- FALSE
liveSampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
invisible(liveSampler$run(30L, 5L))

trainX <- keptSampler$data@x

# routing the training predictors reproduces the default n on every path
expect_equal(keptSampler$getTrees(newdata = trainX), keptSampler$getTrees())
expect_equal(
  keptSampler$getTrees(current = TRUE, newdata = trainX),
  keptSampler$getTrees(current = TRUE)
)
expect_equal(liveSampler$getTrees(newdata = trainX), liveSampler$getTrees())

# independent descent of a flattened tree, counting observations per node
countByDescent <- function(tree, x) {
  counts <- integer(nrow(tree))
  recurse <- function(rows, indices) {
    pos <- rows[1L]
    counts[pos] <<- length(indices)
    if (tree$var[pos] == -1L) {
      return(1L)
    }
    goesLeft <- x[indices, tree$var[pos]] <= tree$value[pos]
    nLeft <- recurse(rows[-1L], indices[goesLeft])
    nRight <- recurse(
      rows[seq.int(2L + nLeft, length(rows))],
      indices[!goesLeft]
    )
    1L + nLeft + nRight
  }
  recurse(seq_len(nrow(tree)), seq_len(nrow(x)))
  counts
}

set.seed(21L)
newX <- matrix(rnorm(40L), ncol = 1L)

# live trees: every observation lands in exactly one leaf, counts match descent
replayed <- liveSampler$getTrees(newdata = newX)
leafSums <- with(replayed[replayed$var == -1L, ], tapply(n, tree, sum))
expect_true(all(leafSums == nrow(newX)))
for (t in unique(replayed$tree)) {
  sub <- replayed[replayed$tree == t, ]
  expect_equal(sub$n, countByDescent(sub, newX))
}

# saved trees route the same way
replayedSaved <- keptSampler$getTrees(newdata = newX, sampleNums = 1L)
for (t in unique(replayedSaved$tree)) {
  sub <- replayedSaved[replayedSaved$tree == t, ]
  expect_equal(sub$n, countByDescent(sub, newX))
}

# a mismatched column count is rejected
expect_error(
  liveSampler$getTrees(newdata = matrix(rnorm(20L), ncol = 2L)),
  "number of columns in 'test' must be equal to that of 'x'"
)

# out-of-range chain/sample/tree indices are rejected, not silently dropped
expect_error(liveSampler$getTrees(chainNums = 999L), "'chainNums' must be in")
expect_error(
  keptSampler$getTrees(sampleNums = 999L),
  "'sampleNums' must be in"
)
expect_error(liveSampler$getTrees(treeNums = 999L), "'treeNums' must be in")

# a fractional chain/sample/tree index is refused, naming the argument,
# rather than silently truncated (coerceOrError's integer branch)
expect_error(
  liveSampler$getTrees(chainNums = 1.5),
  "'chainNums' must be a whole number; got '1.5'",
  fixed = TRUE
)
expect_error(
  liveSampler$getTrees(treeNums = 1.5),
  "'treeNums' must be a whole number; got '1.5'",
  fixed = TRUE
)
expect_error(
  keptSampler$getTrees(sampleNums = 1.5),
  "'sampleNums' must be a whole number; got '1.5'",
  fixed = TRUE
)

# exact ties: an observation equal to a split value routes left (x <= split), a
# boundary random continuous newdata never exercises. Take a tree whose root
# splits and feed values straddling and equal to that split.
liveTrees <- liveSampler$getTrees()
roots <- liveTrees[!duplicated(liveTrees$tree), ]
splitTree <- roots$tree[roots$var != -1L][1L]
expect_false(is.na(splitTree))
s <- liveTrees$value[liveTrees$tree == splitTree][1L]
tieX <- matrix(c(s - 1, s, s, s + 1), ncol = 1L) # 3 of 4 are <= s
tied <- liveSampler$getTrees(newdata = tieX)
tieBlock <- tied[tied$tree == splitTree, ]
expect_equal(tieBlock$n[1L], nrow(tieX)) # root sees all observations
expect_equal(tieBlock$n[2L], sum(tieX[, 1L] <= s)) # left child gets the ties

rm(keptSampler, liveSampler, ctrl, trainX, countByDescent, newX)
rm(replayed, leafSums, replayedSaved, sub, t)
rm(liveTrees, roots, splitTree, s, tieX, tied, tieBlock)
rm(x, y, n)


## ---------------------------------------------------------------------------
## getTrees(forest = ) on a single-forest sampler: NULL (the default), an
## explicit 1L, and a length-1 vector are all bitwise the same read, and an
## invalid forest is refused with the same message resolveForestIndex gives
## getForestFits/getForestAmplitudes/getForestVariableCounts/getLeafPrior.
## ---------------------------------------------------------------------------
set.seed(13L)
n <- 40L
x <- rnorm(n)
y <- x + rnorm(n)
forestSampler <- dbarts::dbarts(
  y ~ x,
  data.frame(x = x, y = y),
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.trees = 5L,
    n.samples = 3L,
    n.burn = 0L,
    updateState = TRUE,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 5L
  )
)
invisible(forestSampler$run(10L, 3L))

defaultForest <- forestSampler$getTrees()
expect_equal(forestSampler$getTrees(forest = NULL), defaultForest)
expect_equal(forestSampler$getTrees(forest = 1L), defaultForest)
expect_equal(forestSampler$getTrees(forest = 1:1), defaultForest)
expect_false("forest" %in% names(defaultForest))
expect_false("chain" %in% names(defaultForest))

expect_error(
  forestSampler$getTrees(forest = 0L),
  "'forest' must be a single positive integer",
  fixed = TRUE
)
expect_error(
  forestSampler$getTrees(forest = -1L),
  "'forest' must be a single positive integer",
  fixed = TRUE
)
expect_error(
  forestSampler$getTrees(forest = NA_integer_),
  "'forest' must be a single positive integer",
  fixed = TRUE
)
expect_error(
  forestSampler$getTrees(forest = 1.5),
  "'forest' must be a whole number",
  fixed = TRUE
)
expect_error(forestSampler$getTrees(forest = 2L), "out of range")

rm(forestSampler, defaultForest, x, y, n)


## ---------------------------------------------------------------------------
## Each forest's OWN tree count, not forest 1's (control@n.trees), governs
## getTrees' default and validation for that forest: a forest() formula term
## defaults to 50 trees regardless of the fit's own n.trees, and a forests =
## declaration can give either forest more trees than the other.
## ---------------------------------------------------------------------------
set.seed(17L)
n <- 80L
a <- rnorm(n)
b <- rnorm(n)
z <- rbinom(n, 1L, 0.5)
y <- a + b + z * (a - b) + rnorm(n, 0, 0.3)
asymDf <- data.frame(a = a, b = b, z = z, y = y)

# bart()'s own default n.trees (75) for forest 1, forest()'s own default (50)
# for forest 2, neither overridden here
asymFit <- dbarts::bart(
  y ~ a + b + z:forest(a + b),
  asymDf,
  keepTrees = TRUE,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 3L,
  n.burn = 3L,
  verbose = FALSE
)
asymTrees <- dbarts::extract(asymFit, "trees")
expect_equal(range(asymTrees$tree[asymTrees$forest == 1L]), c(1L, 75L))
expect_equal(range(asymTrees$tree[asymTrees$forest == 2L]), c(1L, 50L))
# a small n.trees on forest 1 does not truncate forest 2's own, larger count
asymFitSmall <- dbarts::bart(
  y ~ a + b + z:forest(a + b),
  asymDf,
  keepTrees = TRUE,
  n.trees = 7L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 3L,
  n.burn = 3L,
  verbose = FALSE
)
asymSmallForest2 <- dbarts::extract(asymFitSmall, "trees", forest = 2L)
expect_equal(sort(unique(asymSmallForest2$tree)), seq_len(50L))
# plotTree reaches a forest-2 tree index past forest 1's own count
pdf(NULL)
expect_silent(
  plotTree(asymFit, forest = 2L, treeNum = 50L, chainNum = 1L, sampleNum = 1L)
)
dev.off()
# an explicit treeNums is checked against the NAMED forest's own count
expect_error(
  asymFit$fit$getTrees(forest = 2L, treeNums = 1:75),
  "'treeNums' must be in [1, 50] for forest 2",
  fixed = TRUE
)
rm(asymDf, asymFit, asymTrees, asymFitSmall, asymSmallForest2, a, b, z, y, n)

# a forests = sampler whose SECOND forest carries MORE trees than the first
set.seed(19L)
n <- 70L
xMore <- matrix(runif(n * 2L), n, 2L)
zMore <- rbinom(n, 1L, 0.5)
yMore <- xMore[, 1L] + zMore * xMore[, 2L] + rnorm(n, 0, 0.2)
moreControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 8L,
  n.samples = 3L,
  n.burn = 2L,
  keepTrees = TRUE,
  updateState = FALSE,
  verbose = FALSE
)
moreSampler <- dbarts::dbarts(
  xMore,
  yMore,
  forests = list(forest(), forest(basis = ~ factor(zMore), n.trees = 30L)),
  control = moreControl
)
invisible(moreSampler$run(2L, 3L))
moreTrees <- moreSampler$getTrees()
expect_equal(range(moreTrees$tree[moreTrees$forest == 1L]), c(1L, 8L))
expect_equal(range(moreTrees$tree[moreTrees$forest == 2L]), c(1L, 30L))
rm(
  xMore,
  zMore,
  yMore,
  moreControl,
  moreSampler,
  moreTrees,
  n
)


rm(df, testData, treesArgReason)


## ---------------------------------------------------------------------------
## The forest column appears only where there are several forests or a forests
## = declaration; extract(type = "trees") always carries chain, the sampler's
## getTrees only at more than one chain.
## ---------------------------------------------------------------------------
set.seed(31L)
n <- 60L
x1 <- rnorm(n)
z <- factor(rbinom(n, 1L, 0.5))
y <- x1 + rnorm(n)
yBin <- as.integer(x1 + rnorm(n) > 0)
treeDf <- data.frame(y = y, yBin = yBin, x1 = x1, z = z)
treeCols <- c("sample", "tree", "n", "var", "value")

for (nChains in 1:2) {
  chainCols <- if (nChains > 1L) "chain" else NULL
  fits <- list(
    gaussian = dbarts::bart(
      y ~ x1,
      treeDf,
      n.trees = 3L,
      n.samples = 4L,
      n.burn = 2L,
      n.chains = nChains,
      n.threads = 1L,
      keepTrees = TRUE,
      verbose = FALSE
    ),
    probit = dbarts::bart(
      yBin ~ x1,
      treeDf,
      family = "probit",
      n.trees = 3L,
      n.samples = 4L,
      n.burn = 2L,
      n.chains = nChains,
      n.threads = 1L,
      keepTrees = TRUE,
      verbose = FALSE
    ),
    bartBT = dbarts::bartBT(
      x1,
      y,
      ntree = 3L,
      ndpost = 4L,
      nskip = 2L,
      nchain = nChains,
      nthread = 1L,
      keeptrees = TRUE,
      verbose = FALSE
    )
  )
  for (name in names(fits)) {
    fit <- fits[[name]]
    extracted <- dbarts::extract(fit, type = "trees")
    expect_equal(
      colnames(extracted),
      c("chain", treeCols),
      info = paste(name, nChains)
    )
    expect_equal(
      sort(unique(extracted$chain)),
      seq_len(nChains),
      info = paste(name, nChains)
    )
    expect_equal(
      colnames(fit$fit$getTrees()),
      c(chainCols, treeCols),
      info = paste(name, nChains)
    )
    expect_equal(
      colnames(dbarts::extract(fit, type = "trees", current = TRUE)),
      c("chain", "tree", "n", "var", "value"),
      info = paste(name, nChains)
    )
    expect_equal(
      dbarts::extract(fit, type = "trees", forest = 1L),
      extracted,
      info = paste(name, nChains)
    )
    expect_equal(
      fit$fit$getTrees(forest = 1L),
      fit$fit$getTrees(),
      info = paste(name, nChains)
    )
    expect_error(
      fit$fit$getTrees(forest = 2L),
      "out of range",
      info = paste(name, nChains)
    )
  }
}

treeControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 3L,
  n.samples = 3L,
  n.burn = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 3L
)
# a plotTree on a single-forest sampler still needs no forest
pdf(NULL)
expect_silent(fits$gaussian$fit$plotTree(1L, chainNum = 1L, sampleNum = 1L))
dev.off()

declared <- dbarts::dbarts(
  y ~ x1,
  treeDf,
  forests = list(forest()),
  control = treeControl
)
invisible(declared$run())
expect_equal(colnames(declared$getTrees())[1L], "forest")
expect_true(all(declared$getTrees()$forest == 1L))

twoForests <- dbarts::dbarts(
  y ~ x1,
  treeDf,
  forests = list(forest(), forest(basis = ~z)),
  control = treeControl
)
invisible(twoForests$run())
expect_equal(colnames(twoForests$getTrees()), c("forest", treeCols))
expect_equal(sort(unique(twoForests$getTrees()$forest)), 1:2)
expect_error(twoForests$plotTree(1L, chainNum = 1L), "forest required")

# the dbartsSpec() route, which builds the sampler by hand, marks a one-forest
# declaration the same way
specDeclared <- dbarts::dbartsSpec(
  dbarts::dbartsData(y ~ x1, treeDf),
  treeControl,
  forests = list(forest())
)
specSampler <- new(
  "dbartsSampler",
  specDeclared$control,
  specDeclared$model,
  specDeclared$data
)
invisible(specSampler$run())
expect_equal(colnames(specSampler$getTrees()), c("forest", treeCols))

# a one-forest declaration survives save/load and copy()
declaredFile <- tempfile(fileext = ".rds")
saveRDS(declared, declaredFile)
reloaded <- readRDS(declaredFile)
unlink(declaredFile)
expect_equal(colnames(reloaded$getTrees()), c("forest", treeCols))
expect_equal(colnames(declared$copy()$getTrees()), c("forest", treeCols))

# a control taken from another sampler does not carry that sampler's
# configuration: not a forests declaration, and not a variance forest
plainFromDeclared <- dbarts::dbarts(
  y ~ x1,
  treeDf,
  control = declared$control
)
invisible(plainFromDeclared$run())
expect_equal(colnames(plainFromDeclared$getTrees()), treeCols)
expect_null(attr(plainFromDeclared$control, "bartcore.forestsDeclared"))

plain <- dbarts::dbarts(y ~ x1, treeDf, control = treeControl)
plain$setControl(declared$control)
invisible(plain$run())
expect_equal(colnames(plain$getTrees()), treeCols)

varianceDonor <- dbarts::dbarts(
  y ~ x1,
  treeDf,
  variance = varianceForest(n.trees = 5L),
  control = treeControl
)
expect_false(is.null(attr(varianceDonor$control, "bartcore.variance")))
noVariance <- dbarts::dbarts(y ~ x1, treeDf, control = varianceDonor$control)
expect_null(attr(noVariance$control, "bartcore.variance"))
invisible(noVariance$run())
expect_null(noVariance$run()$variance)
plain$setControl(varianceDonor$control)
expect_null(attr(plain$control, "bartcore.variance"))

# a variance forest's trees are not read by getTrees, so no forest column
varianceFit <- dbarts::bart(
  y ~ x1,
  treeDf,
  n.trees = 3L,
  n.samples = 3L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  variance = varianceForest(n.trees = 5L)
)
expect_equal(
  colnames(dbarts::extract(varianceFit, type = "trees")),
  c("chain", treeCols)
)
expect_equal(colnames(varianceFit$fit$getTrees()), treeCols)
