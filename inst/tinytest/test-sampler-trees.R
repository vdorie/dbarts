source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

df <- with(testData, data.frame(x, y))

# test that base bart extracts trees correctly
n.trees <- 3L
n.samples <- 4L
fit <- dbarts::bart(
  y ~ ., df, nthread = 1L, ntree = n.trees, nskip = 0L,
  ndpost = n.samples, keeptrees = TRUE, verbose = FALSE
)
allTrees <- dbarts::extract(fit, "trees")

expect_true(all(c("sample", "tree") %in% colnames(allTrees)))
expect_true(!("chain" %in% colnames(allTrees)))

combinations <- data.frame(
  sample = rep(seq_len(n.samples), each = n.trees),
  tree   = rep(seq_len(n.trees), times = n.samples)
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
  y ~ ., df, group.by = g,
  n.threads = 1L, n.trees = n.trees, n.burn = 0L, n.thin = 1L,
  n.chains = n.chains,
  n.samples = n.samples,
  keepTrees = TRUE, verbose = FALSE
)
allTrees <- dbarts::extract(fit, "trees")

expect_true(all(c("sample", "chain", "tree") %in% colnames(allTrees)))

combinations <- data.frame(
  chain  = rep(seq_len(n.chains), each = n.trees * n.samples),
  sample = rep(rep(seq_len(n.samples), each = n.trees), times = n.chains),
  tree   = rep(rep(seq_len(n.trees), times = n.samples), times = n.chains)
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
inst <- keptSampler$setPredictor(x + rnorm(n, sd = 0.5), "x", forceUpdate = "partial")
liveTrees <- keptSampler$getTrees(current = TRUE)
expect_false(any(liveTrees$var == -1L & liveTrees$n == 0L))
expect_true(length(inst) == n)

rm(keptSampler, liveSampler, makeSampler, saved, current, liveTrees, inst)
rm(x, y, n)


rm(df, testData)

