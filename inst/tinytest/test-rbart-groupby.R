source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart finds group.by
df <- as.data.frame(testData$x)
colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
df$y <- testData$y
df$g <- testData$g
expect_inherits(
  dbarts::rbart_vi(
    y ~ . - g, df, group.by = g,
    n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
    n.trees = 25L, n.threads = 1L, verbose = FALSE
  ),
  "rbart"
)

g <- df$g
df$g <- NULL
expect_inherits(
  dbarts::rbart_vi(
    y ~ . , df, group.by = g,
    n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
    n.trees = 25L, n.threads = 1L, verbose = FALSE
  ),
  "rbart"
)

y <- testData$y
x <- testData$x
expect_inherits(
  dbarts::rbart_vi(
    y ~ x, group.by = g,
    n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
    n.trees = 25L, n.threads = 1L, verbose = FALSE
  ),
  "rbart"
)

rm(x, y, g, df)


# test that works with missing levels
n.train <- 80L
x <- testData$x[seq_len(n.train),]
y <- testData$y[seq_len(n.train)]
g <- factor(testData$g[seq_len(n.train)])

x.test <- testData$x[seq.int(n.train + 1L, nrow(testData$x)),]
g.test <- factor(testData$g[seq.int(n.train + 1L, nrow(testData$x))], levels(g))
levels(g.test)[5L] <- "6"

# check that predict works when we've fit with missing levels
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x, group.by = g, test = x.test, group.by.test = g.test,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean),
  fitted(rbartFit, sample = "test")
)

# check that predicts works for completely new levels
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(
  predict(rbartFit, x.test, g.test, type = "ranef", combineChains = FALSE)
)
expect_equal(
  ranef.pred[,,as.character(1L:4L)],
  rbartFit$ranef[,,as.character(1L:4L)]
)
expect_true(
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[,,5L:26L], c(1L, 2L), sd))
  ) > 0.90
)

# check again with combineChains as TRUE at the top level
g.test <- droplevels(g.test)
levels(g.test) <- c(levels(g)[-5L], "6")
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x, group.by = g, test = x.test, group.by.test = g.test,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = TRUE,
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean),
  fitted(rbartFit, sample = "test")
)

levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(rbartFit, x.test, g.test, type = "ranef"))
expect_equal(
  as.numeric(ranef.pred[,as.character(1L:4L)]),
  as.numeric(rbartFit$ranef[,as.character(1L:4L)])
)
expect_true(
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[,5L:26L], 1L, sd))
  ) > 0.90
)

# check one last time with one chain
g.test <- droplevels(g.test)
levels(g.test) <- c(levels(g)[-5L], "6")
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x, group.by = g, test = x.test, group.by.test = g.test,
  n.samples = 14L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(rbartFit, x.test, g.test, type = "ranef"))
expect_equal(
  as.numeric(ranef.pred[,as.character(1L:4L)]),
  as.numeric(rbartFit$ranef[,as.character(1L:4L)])
)
expect_true(
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[,5L:26L], 1L, sd))
  ) > 0.90
)

# check with more than one missing level
levels(g.test)[4L] <- "7"
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x, group.by = g, test = x.test, group.by.test = g.test,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 4L,
  n.trees = 25L, n.threads = 1L,
  verbose = FALSE
))
expect_inherits(rbartFit, "rbart")

rm(rbartFit, ranef.pred, g.test)

rm(x.test, g, y, x, n.train)

rm(testData)

