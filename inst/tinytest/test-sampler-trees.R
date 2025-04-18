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

rm(df, testData)

