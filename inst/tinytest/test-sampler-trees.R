source(
  system.file("common", "friedmanData.R", package = "dbarts"),
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

expect_true(all(c("sample", "tree") %in% colnames(allTrees)))
expect_true(!("chain" %in% colnames(allTrees)))

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

rm(individualSamples, combinations, allTrees, fit, n.samples, n.trees)


n.g <- 5L
g <- sample(n.g, length(testData$y), replace = TRUE)

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

df$y <- df$y + b[g]
df$g <- g
rm(g, b, sigma.b, n.g)

# test that rbart extracts trees correctly
n.trees <- 3L
n.samples <- 4L
n.chains <- 2L
fit <- dbarts::rbart_vi(
  y ~ .,
  df,
  group.by = g,
  n.threads = 1L,
  n.trees = n.trees,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = n.chains,
  n.samples = n.samples,
  keepTrees = TRUE,
  verbose = FALSE
)
allTrees <- dbarts::extract(fit, "trees")

expect_true(all(c("sample", "chain", "tree") %in% colnames(allTrees)))

combinations <- data.frame(
  chain = rep(seq_len(n.chains), each = n.trees * n.samples),
  sample = rep(rep(seq_len(n.samples), each = n.trees), times = n.chains),
  tree = rep(rep(seq_len(n.trees), times = n.samples), times = n.chains)
)
expect_true(all(
  paste0(combinations$sample, ";", combinations$tree) %in%
    paste0(allTrees$sample, ";", allTrees$tree)
))

individualSamples <- lapply(
  seq_len(n.chains),
  function(i) extract(fit, "trees", chainNums = i)
)
individualSamples <- Reduce(rbind, individualSamples)
row.names(individualSamples) <- as.character(seq_len(nrow(individualSamples)))

expect_equal(allTrees, individualSamples)

rm(individualSamples, combinations, allTrees, fit)
rm(n.chains, n.samples, n.trees)


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
makeSampler <- function(keepTrees) {
  ctrl <- dbarts::dbartsControl(
    n.chains = 1L,
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = TRUE,
    keepTrees = keepTrees,
    verbose = FALSE,
    rngSeed = 9L
  )
  sampler <- dbarts::dbarts(y ~ x, data.frame(x = x, y = y), control = ctrl)
  invisible(sampler$run(30L, 5L))
  sampler
}
keptSampler <- makeSampler(TRUE)
liveSampler <- makeSampler(FALSE)

saved <- keptSampler$getTrees() # snapshots
current <- keptSampler$getTrees(current = TRUE) # live working trees

expect_true("sample" %in% names(saved)) # snapshots carry a sample dimension
expect_false("sample" %in% names(current)) # live trees do not
# the live trees are exactly what an otherwise-identical keepTrees = FALSE
# sampler reports
expect_equal(current, liveSampler$getTrees())
# current = TRUE is ignored (harmlessly) when there are no saved trees
expect_equal(liveSampler$getTrees(current = TRUE), liveSampler$getTrees())

# a partial update keeps the live trees valid; getTrees(current = TRUE) sees it
inst <- keptSampler$setPredictor(
  x + rnorm(n, sd = 0.5),
  "x",
  forceUpdate = "partial"
)
liveTrees <- keptSampler$getTrees(current = TRUE)
expect_false(any(liveTrees$var == -1L & liveTrees$n == 0L))
expect_true(length(inst) == n)

rm(keptSampler, liveSampler, makeSampler, saved, current, liveTrees, inst)
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
  rngSeed = 3L
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
expect_error(liveSampler$getTrees(newdata = matrix(rnorm(20L), ncol = 2L)))

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


rm(df, testData)
